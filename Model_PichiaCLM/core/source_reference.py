"""Dynamic fetch + local cache for source-organism reference data.

Two things are needed to harmonize a candidate's %MinMax profile against a
source organism (see ADR-0005): the source organism's own codon-usage
frequency table, and the *actual* native coding sequence of the specific
gene -- the amino acid sequence alone cannot tell us which synonymous codon
nature used at each position, and that choice is exactly what a %MinMax
harmonization target is built from.

Neither is bundled with the repository or hardcoded: both are fetched from a
public database on first use and cached locally under
``Training/ExternalReferenceCache`` (gitignored -- this is fetched fact, not
a curated training asset). Nothing here silently falls back to host data or
an empty table on failure; every failure path raises so a missing target
never gets confused with a validated one.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Callable

from .analysis import CODON_TO_AA, build_fraction_table
from .biology import normalize_dna

DEFAULT_CACHE_DIR = Path(__file__).resolve().parents[1] / "Training" / "ExternalReferenceCache"
KAZUSA_URL_TEMPLATE = "https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species={taxon_id}&aa=1&style=N"
NCBI_EFETCH_URL_TEMPLATE = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={accession}&rettype=fasta_cds_na&retmode=text"
)

FetchFn = Callable[[str, float], str]


def _default_fetch(url: str, timeout: float) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as response:  # noqa: S310 - fixed https endpoints above
        return response.read().decode("utf-8", errors="replace")


def _parse_kazusa_table(raw_text: str) -> dict[str, int]:
    """Parse Kazusa's ``style=N`` plain-text codon usage table into raw counts.

    Each data line looks like ``UUU F 0.46 17.6(   714298)``: RNA codon,
    amino acid, fraction-within-family, per-thousand, then the absolute
    count in parentheses -- the count is what we need, since fractions are
    recomputed locally with `build_fraction_table` for consistency with how
    the training-data reference is built elsewhere in this module.
    """
    counts: dict[str, int] = {}
    for line in raw_text.splitlines():
        stripped = line.strip()
        if not stripped or "(" not in stripped:
            continue
        before_paren, _, remainder = stripped.partition("(")
        count_text = remainder.split(")")[0].strip()
        fields = before_paren.split()
        if len(fields) < 2 or not count_text.isdigit():
            continue
        rna_codon = fields[0]
        dna_codon = rna_codon.upper().replace("U", "T")
        if dna_codon not in CODON_TO_AA:
            continue
        counts[dna_codon] = int(count_text)
    return counts


def load_source_organism_codon_fractions(
    taxon_id: int,
    *,
    cache_dir: Path | str | None = None,
    fetch: FetchFn = _default_fetch,
    timeout: float = 10.0,
) -> tuple[dict[str, float], int]:
    """Return (fractions, total_codon_count) for a source organism's codon usage.

    Resolution order: local cache for ``taxon_id`` -> fetch from the Kazusa
    Codon Usage Database (same source and query style already used for this
    module's bundled host-side ``PUBLIC_PICHIA_PASTORIS_FRACTIONS``) -> raise.
    Never falls back to host-organism data or an empty table -- see ADR-0005.
    """
    cache_path = _cache_path(cache_dir, f"kazusa_taxon_{taxon_id}.json")
    cached = _read_cache(cache_path)
    if cached is not None:
        return cached["fractions"], cached["total_codon_count"]

    url = KAZUSA_URL_TEMPLATE.format(taxon_id=taxon_id)
    try:
        raw_text = fetch(url, timeout)
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        raise RuntimeError(
            f"Could not fetch codon usage for taxon {taxon_id} from Kazusa and no local "
            f"cache exists at {cache_path}. Not falling back to host-organism data -- "
            "that would silently misrepresent it as a validated source-organism table."
        ) from error

    counts = _parse_kazusa_table(raw_text)
    if not counts:
        raise RuntimeError(f"Kazusa response for taxon {taxon_id} did not contain a parseable codon usage table.")

    fractions = build_fraction_table(Counter(counts))
    total_codon_count = sum(counts.values())
    _write_cache(cache_path, {"fractions": fractions, "total_codon_count": total_codon_count})
    return fractions, total_codon_count


def load_native_source_cds(
    *,
    accession: str | None = None,
    manual_cds: str | None = None,
    cache_dir: Path | str | None = None,
    fetch: FetchFn = _default_fetch,
    timeout: float = 10.0,
) -> str:
    """Return the native nucleotide CDS for a source-organism gene.

    ``manual_cds`` always wins when given -- a researcher who already knows
    exactly which isoform/construct they mean can hand it over directly, no
    network or cache involved (see ADR-0005: SPP1/osteopontin alone has
    multiple transcript variants, so an auto-fetched accession is a default,
    not a guarantee it is the one actually in use). Otherwise resolves via
    local cache for ``accession`` -> fetch from NCBI RefSeq -> raise.
    """
    if manual_cds is not None:
        normalized = normalize_dna(manual_cds)
        if not normalized or any(base not in {"A", "T", "G", "C"} for base in normalized):
            raise ValueError("manual_cds must be a non-empty DNA/RNA sequence using only A/T/G/C/U bases.")
        return normalized

    if not accession:
        raise ValueError("Either manual_cds or accession must be supplied.")

    cache_path = _cache_path(cache_dir, f"{accession}_cds.json")
    cached = _read_cache(cache_path)
    if cached is not None:
        return cached["cds"]

    url = NCBI_EFETCH_URL_TEMPLATE.format(accession=accession)
    try:
        raw_text = fetch(url, timeout)
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        raise RuntimeError(
            f"Could not fetch the native CDS for accession {accession!r} from NCBI and no "
            f"local cache exists at {cache_path}. Pass manual_cds explicitly instead of "
            "relying on a substitute sequence."
        ) from error

    cds = _parse_fasta_cds(raw_text)
    if not cds:
        raise RuntimeError(f"NCBI response for accession {accession!r} did not contain a parseable CDS.")

    _write_cache(cache_path, {"cds": cds})
    return cds


def build_harmonization_target(
    *,
    source_taxon_id: int | None,
    source_native_cds: str | None,
    accession: str | None = None,
    cache_dir: Path | str | None = None,
    fetch: FetchFn = _default_fetch,
    timeout: float = 10.0,
):
    """Assemble a ``MinMaxHarmonizationTarget``, or ``None`` if not requested.

    Lives in ``core`` rather than in each interface so CLI/API/Streamlit share
    one definition of what a valid harmonization request is instead of three.

    Supplying only one half is an error, not a silent no-op: a caller who set
    a taxon id but no native CDS (or vice versa) asked for harmonization and
    must not be handed an un-harmonized ranking that looks like it worked.
    """
    from .candidates import MinMaxHarmonizationTarget

    has_cds_request = bool(source_native_cds) or bool(accession)
    if source_taxon_id is None and not has_cds_request:
        return None
    if source_taxon_id is None:
        raise ValueError(
            "A source native CDS was supplied without source_taxon_id. Harmonization needs both: the CDS "
            "says which synonymous codon the source organism used, the taxon id resolves that organism's "
            "usage table."
        )
    if not has_cds_request:
        raise ValueError(
            "source_taxon_id was supplied without a source native CDS (or accession). The amino acid "
            "sequence cannot substitute -- it does not carry which synonymous codon was used at each "
            "position, which is exactly what harmonization matches against."
        )

    fractions, _ = load_source_organism_codon_fractions(
        source_taxon_id, cache_dir=cache_dir, fetch=fetch, timeout=timeout
    )
    cds = load_native_source_cds(
        accession=accession,
        manual_cds=source_native_cds,
        cache_dir=cache_dir,
        fetch=fetch,
        timeout=timeout,
    )
    return MinMaxHarmonizationTarget(source_cds=cds, source_fractions=fractions)


def _parse_fasta_cds(raw_text: str) -> str:
    lines = [line.strip() for line in raw_text.splitlines()]
    sequence_lines = [line for line in lines if line and not line.startswith(">")]
    return normalize_dna("".join(sequence_lines))


def _cache_path(cache_dir: Path | str | None, file_name: str) -> Path:
    directory = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR
    return directory / file_name


def _read_cache(cache_path: Path) -> dict | None:
    if not cache_path.exists():
        return None
    with cache_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _write_cache(cache_path: Path, payload: dict) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", encoding="utf-8") as handle:
        json.dump({**payload, "fetched_at": time.time()}, handle, ensure_ascii=False, indent=2)

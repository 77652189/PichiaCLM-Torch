from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .analysis import PUBLIC_PICHIA_PASTORIS_FRACTIONS, load_training_codon_reference
from .biology import (
    AA_TO_CODONS,
    CODON_TO_AA,
    DNA_BASES,
    check_translation,
    normalize_amino_acid_sequence,
    normalize_dna,
    split_codons,
    translate_cds,
)


AA_NAME_TO_CODE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "STOP": "*",
}


@dataclass(frozen=True)
class CodonCell:
    codon_number: int
    start: int
    end: int
    dna_codon: str
    rna_codon: str
    amino_acid: str
    synonymous_codons: list[str]
    replaceable: bool
    training_fraction: float | None
    public_fraction: float | None


@dataclass(frozen=True)
class CodonReplacement:
    codon_number: int
    start: int
    end: int
    amino_acid: str
    old_codon: str
    new_codon: str


@dataclass(frozen=True)
class CodonEditResult:
    original_cds: str
    edited_cds: str
    expected_amino_acids: str
    translated_amino_acids: str
    translation_preserved: bool
    replacements: list[CodonReplacement]


def dna_to_rna(codon: str) -> str:
    return normalize_dna(codon).replace("T", "U")


def validate_editable_cds(cds: str) -> str:
    normalized = normalize_dna(cds)
    if not normalized:
        raise ValueError("CDS is empty.")
    invalid_bases = sorted({base for base in normalized if base not in DNA_BASES})
    if invalid_bases:
        raise ValueError(f"CDS contains invalid bases: {', '.join(invalid_bases)}.")
    if len(normalized) % 3 != 0:
        raise ValueError("CDS length must be a multiple of 3.")
    return normalized


def build_codon_cells(cds: str) -> list[CodonCell]:
    normalized = validate_editable_cds(cds)
    codons = split_codons(normalized)
    training_fractions, _ = load_training_codon_reference()
    cells = []
    for index, codon in enumerate(codons, start=1):
        amino_acid = CODON_TO_AA.get(codon, "X")
        synonymous = list(AA_TO_CODONS.get(amino_acid, []))
        replaceable = amino_acid != "*" and len(synonymous) > 1
        cells.append(
            CodonCell(
                codon_number=index,
                start=(index - 1) * 3 + 1,
                end=index * 3,
                dna_codon=codon,
                rna_codon=dna_to_rna(codon),
                amino_acid=amino_acid,
                synonymous_codons=synonymous,
                replaceable=replaceable,
                training_fraction=round(training_fractions.get(codon, 0.0), 4),
                public_fraction=round(PUBLIC_PICHIA_PASTORIS_FRACTIONS.get(codon, 0.0), 4),
            )
        )
    return cells


def search_codon_cells(cells: list[CodonCell], query: str) -> list[int]:
    clean = query.strip().upper()
    if not clean:
        return []

    normalized_dna = normalize_dna(clean)
    if len(normalized_dna) == 3 and all(base in DNA_BASES for base in normalized_dna):
        return [
            cell.codon_number
            for cell in cells
            if cell.dna_codon == normalized_dna
        ]

    amino_acid = normalize_amino_acid_query(clean)
    if amino_acid is None:
        return []
    return [
        cell.codon_number
        for cell in cells
        if cell.amino_acid == amino_acid
    ]


def normalize_amino_acid_query(query: str) -> str | None:
    clean = query.strip().upper()
    if clean in CODON_TO_AA.values():
        return clean
    return AA_NAME_TO_CODE.get(clean)


def codon_options_for_amino_acid(amino_acid: str) -> list[dict[str, object]]:
    aa = normalize_amino_acid_query(amino_acid)
    if aa is None:
        return []
    training_fractions, _ = load_training_codon_reference()
    rows = []
    for codon in AA_TO_CODONS.get(aa, []):
        rows.append(
            {
                "codon": codon,
                "rna_codon": dna_to_rna(codon),
                "amino_acid": aa,
                "training_fraction": round(training_fractions.get(codon, 0.0), 4),
                "public_fraction": round(PUBLIC_PICHIA_PASTORIS_FRACTIONS.get(codon, 0.0), 4),
            }
        )
    return rows


def replace_selected_codons(
    cds: str,
    selected_indices: Iterable[int],
    target_codon: str,
    expected_amino_acids: str | None = None,
) -> CodonEditResult:
    normalized = validate_editable_cds(cds)
    codons = split_codons(normalized)
    selected = sorted({int(index) for index in selected_indices})
    if not selected:
        raise ValueError("Select at least one codon to replace.")
    if selected[0] < 1 or selected[-1] > len(codons):
        raise ValueError("Selected codon number is outside the CDS range.")

    selected_amino_acids = {
        CODON_TO_AA.get(codons[index - 1], "X")
        for index in selected
    }
    if len(selected_amino_acids) != 1:
        raise ValueError("Selected codons encode different amino acids; replace one amino acid group at a time.")
    amino_acid = next(iter(selected_amino_acids))
    if amino_acid == "*" or len(AA_TO_CODONS.get(amino_acid, [])) <= 1:
        raise ValueError(f"Amino acid {amino_acid} has no replaceable synonymous codon.")

    target = normalize_dna(target_codon)
    if len(target) != 3 or target not in CODON_TO_AA:
        raise ValueError(f"Target codon is not valid: {target_codon}.")
    if CODON_TO_AA[target] != amino_acid:
        raise ValueError(
            f"Target codon {target} encodes {CODON_TO_AA[target]}, not selected amino acid {amino_acid}."
        )

    edited_codons = codons.copy()
    replacements = []
    for codon_number in selected:
        old_codon = edited_codons[codon_number - 1]
        if old_codon == target:
            continue
        edited_codons[codon_number - 1] = target
        replacements.append(
            CodonReplacement(
                codon_number=codon_number,
                start=(codon_number - 1) * 3 + 1,
                end=codon_number * 3,
                amino_acid=amino_acid,
                old_codon=old_codon,
                new_codon=target,
            )
        )

    edited_cds = "".join(edited_codons)
    expected = (
        normalize_amino_acid_sequence(expected_amino_acids)
        if expected_amino_acids
        else translate_cds(normalized)
    )
    translation_check = check_translation(edited_cds, expected)
    if translation_check.translation_matches_input is not True:
        raise ValueError("Replacement changed the translated amino acid sequence.")

    return CodonEditResult(
        original_cds=normalized,
        edited_cds=edited_cds,
        expected_amino_acids=expected,
        translated_amino_acids=translation_check.translated_amino_acids,
        translation_preserved=True,
        replacements=replacements,
    )

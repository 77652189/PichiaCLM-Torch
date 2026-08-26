from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from Model_PichiaCLM.core.source_reference import (
    build_harmonization_target,
    load_native_source_cds,
    load_source_organism_codon_fractions,
)

FAKE_KAZUSA_TEXT = """
UUU F 0.46 17.6(  714298)
UUC F 0.54 20.3(  824692)
GGG G 0.10   1.0(   1000)
GGA G 0.33   3.0(   3300)
GGT G 0.44   4.0(   4400)
GGC G 0.14   1.5(   1500)
"""


def _raising_fetch(url: str, timeout: float) -> str:
    raise OSError("network should not have been called")


class SourceOrganismCodonFractionsTests(unittest.TestCase):
    def test_fetches_parses_and_caches_when_no_cache_exists(self) -> None:
        with tempfile.TemporaryDirectory() as cache_dir:
            calls = []

            def fake_fetch(url: str, timeout: float) -> str:
                calls.append(url)
                return FAKE_KAZUSA_TEXT

            fractions, total_count = load_source_organism_codon_fractions(
                9606, cache_dir=cache_dir, fetch=fake_fetch
            )

            self.assertEqual(len(calls), 1)
            self.assertIn("species=9606", calls[0])
            self.assertEqual(total_count, 714298 + 824692 + 1000 + 3300 + 4400 + 1500)
            self.assertAlmostEqual(fractions["TTT"], 714298 / (714298 + 824692))
            self.assertAlmostEqual(fractions["GGG"], 1000 / (1000 + 3300 + 4400 + 1500))

            cache_file = Path(cache_dir) / "kazusa_taxon_9606.json"
            self.assertTrue(cache_file.exists())
            cached_payload = json.loads(cache_file.read_text(encoding="utf-8"))
            self.assertEqual(cached_payload["total_codon_count"], total_count)

    def test_uses_cache_without_calling_fetch(self) -> None:
        with tempfile.TemporaryDirectory() as cache_dir:
            cache_file = Path(cache_dir) / "kazusa_taxon_9606.json"
            cache_file.write_text(
                json.dumps({"fractions": {"TTT": 0.5, "TTC": 0.5}, "total_codon_count": 100}),
                encoding="utf-8",
            )

            fractions, total_count = load_source_organism_codon_fractions(
                9606, cache_dir=cache_dir, fetch=_raising_fetch
            )

            self.assertEqual(fractions, {"TTT": 0.5, "TTC": 0.5})
            self.assertEqual(total_count, 100)

    def test_raises_instead_of_falling_back_when_fetch_fails_and_no_cache(self) -> None:
        with tempfile.TemporaryDirectory() as cache_dir:
            with self.assertRaisesRegex(RuntimeError, "Could not fetch codon usage"):
                load_source_organism_codon_fractions(9606, cache_dir=cache_dir, fetch=_raising_fetch)


class NativeSourceCdsTests(unittest.TestCase):
    def test_manual_cds_wins_and_never_touches_the_network(self) -> None:
        result = load_native_source_cds(manual_cds="atg aaa ttt\ngggtaa", fetch=_raising_fetch)
        self.assertEqual(result, "ATGAAATTTGGGTAA")

    def test_manual_cds_rejects_non_dna_characters(self) -> None:
        with self.assertRaises(ValueError):
            load_native_source_cds(manual_cds="ATGXYZ", fetch=_raising_fetch)

    def test_requires_either_manual_cds_or_accession(self) -> None:
        with self.assertRaises(ValueError):
            load_native_source_cds(fetch=_raising_fetch)

    def test_fetches_parses_fasta_and_caches_when_no_cache_exists(self) -> None:
        with tempfile.TemporaryDirectory() as cache_dir:
            calls = []

            def fake_fetch(url: str, timeout: float) -> str:
                calls.append(url)
                return ">NM_999999.1 Homo sapiens fake gene, mRNA\nATGAAATTT\nGGGTAA\n"

            cds = load_native_source_cds(accession="NM_999999", cache_dir=cache_dir, fetch=fake_fetch)

            self.assertEqual(cds, "ATGAAATTTGGGTAA")
            self.assertEqual(len(calls), 1)
            self.assertIn("id=NM_999999", calls[0])
            self.assertTrue((Path(cache_dir) / "NM_999999_cds.json").exists())

    def test_uses_cache_without_calling_fetch(self) -> None:
        with tempfile.TemporaryDirectory() as cache_dir:
            (Path(cache_dir) / "NM_999999_cds.json").write_text(
                json.dumps({"cds": "ATGAAATAA"}), encoding="utf-8"
            )

            cds = load_native_source_cds(accession="NM_999999", cache_dir=cache_dir, fetch=_raising_fetch)

            self.assertEqual(cds, "ATGAAATAA")

    def test_raises_instead_of_falling_back_when_fetch_fails_and_no_cache(self) -> None:
        with tempfile.TemporaryDirectory() as cache_dir:
            with self.assertRaisesRegex(RuntimeError, "Could not fetch the native CDS"):
                load_native_source_cds(accession="NM_999999", cache_dir=cache_dir, fetch=_raising_fetch)


if __name__ == "__main__":
    unittest.main()


class BuildHarmonizationTargetTests(unittest.TestCase):
    """Partial configuration must fail loudly: a half-specified harmonization
    request that silently produced no target would hand back an un-harmonized
    ranking the caller believes is harmonized."""

    def test_returns_none_when_nothing_requested(self) -> None:
        self.assertIsNone(
            build_harmonization_target(source_taxon_id=None, source_native_cds=None)
        )

    def test_cds_without_taxon_id_raises(self) -> None:
        with self.assertRaisesRegex(ValueError, "without source_taxon_id"):
            build_harmonization_target(source_taxon_id=None, source_native_cds="ATGGCT")

    def test_taxon_id_without_cds_raises(self) -> None:
        with self.assertRaisesRegex(ValueError, "without a source native CDS"):
            build_harmonization_target(source_taxon_id=9606, source_native_cds=None)

    def test_builds_target_from_both_halves_without_network(self) -> None:
        def fake_fetch(url: str, timeout: float) -> str:
            assert "9606" in url
            return FAKE_KAZUSA_TEXT

        with tempfile.TemporaryDirectory() as cache_dir:
            target = build_harmonization_target(
                source_taxon_id=9606,
                source_native_cds="ATGGCTGCG",
                cache_dir=cache_dir,
                fetch=fake_fetch,
            )
        self.assertEqual(target.source_cds, "ATGGCTGCG")
        self.assertIn("GGT", target.source_fractions)

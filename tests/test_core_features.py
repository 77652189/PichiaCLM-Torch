from __future__ import annotations

import unittest

import torch

from Model_PichiaCLM.core.analysis import analyze_cds, load_training_codon_reference
from Model_PichiaCLM.core.biology import check_translation, translate_cds
from Model_PichiaCLM.core.candidates import (
    CandidateGenerationOptions,
    codon_preference_stats,
    compare_cds,
    generate_cds_candidates,
)
from Model_PichiaCLM.core.fasta import FastaRecord, format_fasta, parse_fasta
from Model_PichiaCLM.core.fusion import compare_signal_fusion
from Model_PichiaCLM.core.postprocess import conservative_postprocess
from Model_PichiaCLM.core.restriction import scan_restriction_sites
from Model_PichiaCLM.core.schemas import PredictionResult
from Model_PichiaCLM.interfaces.api import AnalyzeCdsRequest, CdsRecord, AnalyzeCdsBatchRequest, analyze_cds_batch, analyze_cds_endpoint


class CoreBiologyTests(unittest.TestCase):
    def test_translate_and_check_known_sequence(self) -> None:
        cds = "ATGTCCACAAATCCCAAACCACAGAGA"
        self.assertEqual(translate_cds(cds), "MSTNPKPQR")
        check = check_translation(cds, "MSTNPKPQR")
        self.assertTrue(check.valid_dna)
        self.assertTrue(check.length_multiple_of_three)
        self.assertTrue(check.translation_matches_input)
        self.assertEqual(check.internal_stop_codons, [])

    def test_translation_check_flags_frame_stop_and_invalid_base(self) -> None:
        invalid = check_translation("ATGTAGGCTN", "MAA")
        self.assertFalse(invalid.valid_dna)
        self.assertIn("N", invalid.invalid_bases)
        self.assertFalse(invalid.length_multiple_of_three)
        self.assertEqual(invalid.internal_stop_codons, [2])
        self.assertFalse(invalid.translation_matches_input)


class FastaTests(unittest.TestCase):
    def test_parse_and_format_batch_fasta(self) -> None:
        records = parse_fasta(">seq1 first\nMSTN\nPK\n>seq2\nMST\n")
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].description, "first")
        self.assertEqual(records[0].sequence, "MSTNPK")
        rendered = format_fasta([FastaRecord("cds1", "optimized", "ATG" * 30)], line_width=9)
        self.assertIn(">cds1 optimized", rendered)
        self.assertIn("ATGATGATG", rendered)


class AnalysisAndRestrictionTests(unittest.TestCase):
    def test_analysis_reports_gc_cai_and_restriction_sites(self) -> None:
        report = analyze_cds("GAATTC", amino_acids="EF")
        self.assertEqual(report.gc_percent, 33.33)
        self.assertIsNotNone(report.cai.training)
        self.assertIsNotNone(report.cai.public)
        self.assertTrue(report.translation_matches_input)
        self.assertEqual(report.restriction_sites[0].name, "EcoRI")

    def test_restriction_scan_supports_custom_sites(self) -> None:
        hits = scan_restriction_sites("AAACCCAAACCC", include_defaults=False, custom_sites={"Custom": "AAACCC"})
        self.assertEqual(len(hits), 2)
        self.assertEqual([hit.start for hit in hits], [1, 7])

    def test_api_analyzes_external_cds_without_prediction(self) -> None:
        payload = analyze_cds_endpoint(
            AnalyzeCdsRequest(
                cds="ATGTCCACAAATCCCAAACCACAGAGA",
                expected_amino_acids="MSTNPKPQR",
            )
        )
        self.assertEqual(payload["translated_amino_acids"], "MSTNPKPQR")
        self.assertTrue(payload["analysis"]["translation_matches_input"])

    def test_api_analyzes_external_cds_batch(self) -> None:
        payload = analyze_cds_batch(
            AnalyzeCdsBatchRequest(
                records=[
                    CdsRecord(
                        id="cds1",
                        cds="ATGTCCACAAATCCCAAACCACAGAGA",
                        expected_amino_acids="MSTNPKPQR",
                    )
                ]
            )
        )
        self.assertEqual(payload["records"][0]["id"], "cds1")
        self.assertTrue(payload["records"][0]["analysis"]["translation_matches_input"])


class PostprocessTests(unittest.TestCase):
    def test_postprocess_removes_ecori_without_changing_translation(self) -> None:
        reference, _ = load_training_codon_reference()
        result = conservative_postprocess("GAATTC", "EF", reference_fractions=reference)
        self.assertTrue(result.translation_preserved)
        self.assertEqual(result.amino_acids, "EF")
        self.assertNotIn("GAATTC", result.optimized_cds)
        self.assertGreaterEqual(len(result.replacements), 1)

    def test_postprocess_attempts_to_reduce_high_local_gc(self) -> None:
        reference, _ = load_training_codon_reference()
        original = "CCG" * 10
        result = conservative_postprocess(original, "P" * 10, reference_fractions=reference)
        self.assertTrue(result.translation_preserved)
        self.assertNotEqual(result.optimized_cds, original)
        self.assertLess(result.optimized_cds.count("G") + result.optimized_cds.count("C"), original.count("G") + original.count("C"))


class FusionTests(unittest.TestCase):
    def test_signal_fusion_compares_whole_and_segmented_modes(self) -> None:
        predictor = FakePredictor()
        comparison = compare_signal_fusion(predictor, signal_peptide="MS", mature_protein="TP")
        self.assertEqual(comparison.fused_amino_acids, "MSTP")
        self.assertTrue(comparison.whole_sequence.analysis.translation_matches_input)
        self.assertTrue(comparison.segmented.analysis.translation_matches_input)
        self.assertEqual(comparison.whole_sequence.cleavage_window.amino_acids, "MSTP")


class CandidateGenerationTests(unittest.TestCase):
    def test_candidates_translate_to_input_and_are_reproducible(self) -> None:
        predictor = FakePredictor()
        options = CandidateGenerationOptions(num_candidates=5, seed=42, max_attempts=100)
        first = generate_cds_candidates(predictor, "MSTNPKPQR", options=options)
        second = generate_cds_candidates(predictor, "MSTNPKPQR", options=options)

        self.assertEqual([candidate.cds for candidate in first.candidates], [candidate.cds for candidate in second.candidates])
        self.assertLessEqual(first.generated_candidates, 5)
        self.assertGreater(len({candidate.cds for candidate in first.candidates}), 1)
        for candidate in first.candidates:
            self.assertEqual(candidate.analysis.translated_amino_acids, "MSTNPKPQR")
            self.assertTrue(candidate.analysis.translation_matches_input)
            self.assertLessEqual(candidate.difference_from_reference.codon_difference_percent, 20.0)
            self.assertLessEqual(
                candidate.codon_preference.avoidable_lowest_count,
                first.candidates[0].codon_preference.avoidable_lowest_count,
            )

    def test_candidate_generation_reports_exhausted_design_space(self) -> None:
        predictor = FakePredictor()
        result = generate_cds_candidates(
            predictor,
            "MWMWM",
            options=CandidateGenerationOptions(num_candidates=3, seed=7, max_attempts=5),
        )

        self.assertEqual(result.generated_candidates, 1)
        self.assertTrue(result.exhausted)
        self.assertIsNotNone(result.note)

    def test_cds_difference_metrics(self) -> None:
        difference = compare_cds("AAACCC", "AAAGGG")
        self.assertEqual(difference.bp_differences, 3)
        self.assertEqual(difference.bp_difference_percent, 50.0)
        self.assertEqual(difference.codon_differences, 1)
        self.assertEqual(difference.codon_difference_percent, 50.0)

    def test_kazusa_codon_preference_stats(self) -> None:
        stats = codon_preference_stats("GCTGCCGCG")
        self.assertEqual(stats.codon_count, 3)
        self.assertEqual(stats.top_preferred_count, 1)
        self.assertEqual(stats.second_preferred_count, 1)
        self.assertEqual(stats.lowest_preferred_count, 1)
        self.assertEqual(stats.avoidable_lowest_count, 1)


class FakePredictor:
    device = torch.device("cpu")
    codons = {
        "M": ["ATG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "N": ["AAT", "AAC"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "K": ["AAA", "AAG"],
        "Q": ["CAA", "CAG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "E": ["GAA", "GAG"],
        "F": ["TTC", "TTT"],
        "W": ["TGG"],
    }

    def predict(self, amino_acids: str, allow_unknown: bool = False) -> PredictionResult:
        cds = "".join(self.codons[aa][0] for aa in amino_acids)
        return PredictionResult(
            amino_acids=amino_acids,
            cds=cds,
            codon_ids=list(range(1, len(amino_acids) + 1)),
            device="test",
        )

    def predict_sample(
        self,
        amino_acids: str,
        allow_unknown: bool = False,
        temperature: float = 0.8,
        generator: torch.Generator | None = None,
    ) -> PredictionResult:
        codons = []
        codon_ids = []
        for position, aa in enumerate(amino_acids, start=1):
            choices = self.codons[aa]
            if len(choices) == 1:
                selected = 0
            else:
                weights = torch.ones(len(choices)) / temperature
                selected = torch.multinomial(weights, 1, generator=generator).item()
            codons.append(choices[selected])
            codon_ids.append(position * 10 + selected)
        return PredictionResult(
            amino_acids=amino_acids,
            cds="".join(codons),
            codon_ids=codon_ids,
            device="test",
        )


if __name__ == "__main__":
    unittest.main()

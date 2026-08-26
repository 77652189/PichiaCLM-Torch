from __future__ import annotations

import unittest
from itertools import combinations

import torch

from Model_PichiaCLM.core.analysis import (
    PUBLIC_PICHIA_PASTORIS_FRACTIONS,
    MIN_MAX_WINDOW,
    analyze_cds,
    compare_min_max_profiles,
    load_training_codon_reference,
    min_max_profile,
)
from Model_PichiaCLM.core.biology import check_translation, split_codons, translate_cds
from Model_PichiaCLM.core.candidates import (
    PLACEHOLDER_MAX_CODON_SIMILARITY_PERCENT,
    STRATEGY_TEMPERATURE_SAMPLING,
    CandidateSubsetSelection,
    MinMaxHarmonizationTarget,
    _rank_subset_by_min_max,
    CandidateGenerationOptions,
    CdsCandidate,
    codon_preference_stats,
    compare_cds,
    generate_cds_candidates,
    pairwise_similarity_rows,
    quality_summary,
    select_low_similarity_subset,
)
from Model_PichiaCLM.core.codon_editor import (
    build_codon_cells,
    replace_selected_codons,
    search_codon_cells,
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


class CodonEditorTests(unittest.TestCase):
    def test_build_codon_cells_tracks_order_and_positions(self) -> None:
        cells = build_codon_cells("ATGTCCACAAATCCCAAACCACAGAGA")

        self.assertEqual(len(cells), 9)
        self.assertEqual(cells[0].codon_number, 1)
        self.assertEqual(cells[0].start, 1)
        self.assertEqual(cells[0].end, 3)
        self.assertEqual(cells[0].dna_codon, "ATG")
        self.assertEqual(cells[0].rna_codon, "AUG")
        self.assertEqual(cells[0].amino_acid, "M")
        self.assertFalse(cells[0].replaceable)
        self.assertEqual(cells[-1].codon_number, 9)
        self.assertEqual(cells[-1].start, 25)
        self.assertEqual(cells[-1].end, 27)
        self.assertEqual(cells[-1].dna_codon, "AGA")

    def test_search_accepts_rna_codon_and_amino_acid_query(self) -> None:
        cells = build_codon_cells("TTCTTTCCACCACCT")

        self.assertEqual(search_codon_cells(cells, "UUC"), [1])
        self.assertEqual(search_codon_cells(cells, "CCA"), [3, 4])
        self.assertEqual(search_codon_cells(cells, "P"), [3, 4, 5])
        self.assertEqual(search_codon_cells(cells, "Phe"), [1, 2])

    def test_replace_selected_codons_preserves_translation_and_analysis(self) -> None:
        result = replace_selected_codons("CCACCACCT", [1, 2], "CCT", expected_amino_acids="PPP")

        self.assertEqual(result.edited_cds, "CCTCCTCCT")
        self.assertTrue(result.translation_preserved)
        self.assertEqual(len(result.replacements), 2)
        report = analyze_cds(result.edited_cds, amino_acids="PPP")
        self.assertTrue(report.translation_matches_input)

    def test_replace_selected_codons_rejects_non_synonymous_target(self) -> None:
        with self.assertRaisesRegex(ValueError, "not selected amino acid"):
            replace_selected_codons("CCACCACCT", [1, 2], "TTC", expected_amino_acids="PPP")

    def test_replace_selected_codons_rejects_mixed_amino_acid_selection(self) -> None:
        with self.assertRaisesRegex(ValueError, "different amino acids"):
            replace_selected_codons("CCATTC", [1, 2], "CCT")

    def test_single_codon_amino_acids_are_not_replaceable(self) -> None:
        cells = build_codon_cells("ATGTGG")

        self.assertFalse(cells[0].replaceable)
        self.assertFalse(cells[1].replaceable)
        with self.assertRaisesRegex(ValueError, "no replaceable"):
            replace_selected_codons("ATGTGG", [1], "ATG")


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
        options = CandidateGenerationOptions(num_candidates=5, seed=42, max_attempts=100, subset_size=3)
        first = generate_cds_candidates(predictor, "MSTNPKPQR", options=options)
        second = generate_cds_candidates(predictor, "MSTNPKPQR", options=options)

        self.assertEqual([candidate.cds for candidate in first.candidates], [candidate.cds for candidate in second.candidates])
        self.assertLessEqual(first.generated_candidates, 5)
        self.assertGreater(len({candidate.cds for candidate in first.candidates}), 1)
        expected_pairs = first.generated_candidates * (first.generated_candidates - 1) // 2
        self.assertEqual(len(first.pairwise_similarities), expected_pairs)
        self.assertEqual(first.pairwise_diversity.comparisons, expected_pairs)
        self.assertIsNotNone(first.pairwise_diversity.mean_codon_difference_percent)
        self.assertIsNotNone(first.recommended_subset)
        self.assertEqual(first.recommended_subset.selected_size, 3)
        self.assertEqual(first.recommended_subset.comparisons, 3)
        candidate_ranks = {candidate.rank for candidate in first.candidates}
        self.assertTrue(set(first.recommended_subset.selected_ranks).issubset(candidate_ranks))
        pair_lookup = {
            tuple(sorted((row.left_rank, row.right_rank))): row
            for row in first.pairwise_similarities
        }
        best_possible_max_similarity = min(
            max(pair_lookup[tuple(sorted(pair))].codon_similarity_percent for pair in combinations(ranks, 2))
            for ranks in combinations(sorted(candidate_ranks), 3)
        )
        self.assertEqual(first.recommended_subset.max_codon_similarity_percent, best_possible_max_similarity)
        for candidate in first.candidates:
            self.assertEqual(candidate.analysis.translated_amino_acids, "MSTNPKPQR")
            self.assertTrue(candidate.analysis.translation_matches_input)
            self.assertLessEqual(candidate.difference_from_reference.codon_difference_percent, 20.0)
            if candidate.rank > 1:
                self.assertEqual(candidate.source, "kazusa_diverse")
                self.assertGreaterEqual(candidate.difference_from_reference.codon_difference_percent, 10.0)
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

    def test_similarity_subset_reports_unmet_constraint_instead_of_silently_passing(self) -> None:
        reference_cds = "CCT" * 6
        candidates = [
            _make_candidate(rank=1, cds=reference_cds, reference_cds=reference_cds),
            _make_candidate(rank=2, cds="CCT" * 5 + "CCC", reference_cds=reference_cds),
        ]
        pairwise = pairwise_similarity_rows(candidates)
        self.assertEqual(pairwise[0].codon_similarity_percent, 83.33)

        strict = select_low_similarity_subset(
            candidates,
            pairwise,
            subset_size=2,
            max_codon_similarity_percent=10.0,
        )
        self.assertFalse(strict.constraint_satisfied)
        self.assertFalse(strict.threshold_is_placeholder)
        self.assertEqual(strict.codon_similarity_threshold_percent, 10.0)
        self.assertEqual(strict.max_codon_similarity_percent, 83.33)

        lenient = select_low_similarity_subset(
            candidates,
            pairwise,
            subset_size=2,
            max_codon_similarity_percent=90.0,
        )
        self.assertTrue(lenient.constraint_satisfied)

    def test_high_bp_similarity_alone_does_not_fail_the_gate(self) -> None:
        """ADR-0007: the gate is codon-axis only.

        This pair is in the same regime as real hLF candidates -- a minority of
        codons changed, each at the wobble base -- which puts bp similarity far
        above any useful ceiling while codon similarity stays low. Gating on bp
        made the constraint permanently unsatisfiable, so it must not come back.
        """
        reference_cds = "CCT" * 8
        candidates = [
            _make_candidate(rank=1, cds=reference_cds, reference_cds=reference_cds),
            _make_candidate(rank=2, cds="CCA" * 2 + "CCT" * 6, reference_cds=reference_cds),
        ]
        pairwise = pairwise_similarity_rows(candidates)
        self.assertEqual(pairwise[0].codon_similarity_percent, 75.0)
        self.assertEqual(pairwise[0].bp_similarity_percent, 91.67)

        subset = select_low_similarity_subset(
            candidates, pairwise, subset_size=2, max_codon_similarity_percent=80.0
        )

        self.assertTrue(
            subset.constraint_satisfied,
            "codon similarity 75% is under the 80% ceiling, so bp similarity 91.67% must not fail it",
        )
        self.assertEqual(subset.max_bp_similarity_percent, 91.67, "bp similarity is still reported")
        self.assertFalse(hasattr(subset, "bp_similarity_threshold_percent"))

    def test_similarity_subset_flags_placeholder_threshold_when_not_supplied(self) -> None:
        reference_cds = "CCT" * 6
        candidates = [
            _make_candidate(rank=1, cds=reference_cds, reference_cds=reference_cds),
            _make_candidate(rank=2, cds="CCT" * 5 + "CCC", reference_cds=reference_cds),
        ]
        pairwise = pairwise_similarity_rows(candidates)

        default = select_low_similarity_subset(candidates, pairwise, subset_size=2)

        self.assertTrue(default.threshold_is_placeholder)
        self.assertEqual(default.codon_similarity_threshold_percent, PLACEHOLDER_MAX_CODON_SIMILARITY_PERCENT)
        self.assertFalse(default.constraint_satisfied)


class MinMaxHarmonizationRankingTests(unittest.TestCase):
    """Ranking is exercised through the module-level helper rather than through
    ``generate_cds_candidates`` because ``FakePredictor`` samples synonymous
    codons uniformly (see the note on its ``predict_sample``), so it trips the
    ADR-0001 avoidable-lowest gate on almost every draw once a sequence is long
    enough to fill a %MinMax window. That is a property of the test double, not
    of the real model -- measurements with the shipped weights are in
    docs/EXECUTION_PLAN.md. Driving the helper directly keeps this test about
    the ranking criterion instead of about the fake's draw distribution.
    """

    def _subset(self, ranks: list[int]) -> CandidateSubsetSelection:
        return CandidateSubsetSelection(
            requested_size=len(ranks),
            selected_size=len(ranks),
            selected_ranks=ranks,
            method="test",
            comparisons=1,
            min_bp_similarity_percent=0.0,
            mean_bp_similarity_percent=0.0,
            max_bp_similarity_percent=0.0,
            min_codon_similarity_percent=0.0,
            mean_codon_similarity_percent=0.0,
            max_codon_similarity_percent=0.0,
            codon_similarity_threshold_percent=80.0,
            threshold_is_placeholder=True,
            constraint_satisfied=True,
        )

    def test_candidate_matching_the_source_shape_is_ranked_first(self) -> None:
        source_cds = "GCG" * 20 + "GCT" * 20
        mirrors_source = "GCG" * 20 + "GCT" * 20
        inverts_source = "GCT" * 20 + "GCG" * 20
        candidates = [
            _make_candidate(rank=1, cds=inverts_source, reference_cds=source_cds),
            _make_candidate(rank=2, cds=mirrors_source, reference_cds=source_cds),
        ]
        target = MinMaxHarmonizationTarget(
            source_cds=source_cds,
            source_fractions=PUBLIC_PICHIA_PASTORIS_FRACTIONS,
        )

        ranked = _rank_subset_by_min_max(self._subset([1, 2]), candidates, harmonization_target=target)

        self.assertEqual(ranked.selected_ranks, [2, 1])
        self.assertEqual(ranked.selected_size, 2, "ranking must reorder only, never drop candidates")

    def test_harmonization_target_is_honored_under_the_default_strategy(self) -> None:
        """A supplied target must not be silently dropped just because the
        default generation strategy is in use -- the caller would get an
        unranked subset that looks ranked."""
        predictor = FakePredictor()
        target = MinMaxHarmonizationTarget(
            source_cds="CCG" * 20,
            source_fractions=PUBLIC_PICHIA_PASTORIS_FRACTIONS,
        )
        result = generate_cds_candidates(
            predictor,
            "MSTNPKPQR",
            options=CandidateGenerationOptions(num_candidates=4, subset_size=3, seed=5),
            harmonization_target=target,
        )
        self.assertEqual(result.recommended_subset.ranking_criterion, "harmonization")

    def test_default_strategy_without_target_reports_no_reranking(self) -> None:
        predictor = FakePredictor()
        result = generate_cds_candidates(
            predictor,
            "MSTNPKPQR",
            options=CandidateGenerationOptions(num_candidates=4, subset_size=3, seed=5),
        )
        self.assertEqual(result.recommended_subset.ranking_criterion, "none")

    def test_temperature_sampling_without_target_uses_host_dip_criterion(self) -> None:
        predictor = FakePredictor()
        result = generate_cds_candidates(
            predictor,
            "MSTNPKPQR",
            options=CandidateGenerationOptions(
                num_candidates=4,
                subset_size=3,
                seed=5,
                strategy=STRATEGY_TEMPERATURE_SAMPLING,
            ),
        )
        self.assertEqual(result.recommended_subset.ranking_criterion, "host_worst_dip")

    def test_without_a_target_it_falls_back_to_the_host_only_dip_criterion(self) -> None:
        source_cds = "GCG" * 20 + "GCT" * 20
        candidates = [
            _make_candidate(rank=1, cds="GCT" * 20 + "GCG" * 20, reference_cds=source_cds),
            _make_candidate(rank=2, cds="GCG" * 20 + "GCT" * 20, reference_cds=source_cds),
        ]

        ranked = _rank_subset_by_min_max(self._subset([1, 2]), candidates)

        self.assertEqual(sorted(ranked.selected_ranks), [1, 2])


class MinMaxProfileComparisonTests(unittest.TestCase):
    def test_identical_profiles_have_zero_distance(self) -> None:
        profile = min_max_profile(["GCT"] * MIN_MAX_WINDOW, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        comparison = compare_min_max_profiles(profile, profile)
        self.assertEqual(comparison.mean_absolute_difference, 0.0)
        self.assertEqual(comparison.max_absolute_difference, 0.0)
        self.assertEqual(comparison.comparable_windows, len(profile))
        self.assertEqual(comparison.skipped_windows, 0)

    def test_opposite_profiles_are_maximally_distant(self) -> None:
        fastest = min_max_profile(["GCT"] * MIN_MAX_WINDOW, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        slowest = min_max_profile(["GCG"] * MIN_MAX_WINDOW, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        comparison = compare_min_max_profiles(fastest, slowest)
        self.assertEqual(comparison.mean_absolute_difference, 200.0)

    def test_length_mismatch_raises_instead_of_aligning_silently(self) -> None:
        short = min_max_profile(["GCT"] * MIN_MAX_WINDOW, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        long = min_max_profile(["GCT"] * (MIN_MAX_WINDOW + 3), PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        with self.assertRaisesRegex(ValueError, "different window counts"):
            compare_min_max_profiles(short, long)

    def test_windows_without_comparable_codons_are_skipped_not_counted_as_agreement(self) -> None:
        untranslatable = min_max_profile(["ATG"] * MIN_MAX_WINDOW, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        real = min_max_profile(["GCT"] * MIN_MAX_WINDOW, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        comparison = compare_min_max_profiles(untranslatable, real)
        self.assertEqual(comparison.comparable_windows, 0)
        self.assertEqual(comparison.skipped_windows, len(real))
        self.assertIsNone(comparison.mean_absolute_difference)


class MinMaxProfileTests(unittest.TestCase):
    def test_all_top_preferred_codons_score_near_positive_100(self) -> None:
        codons = ["GCT"] * MIN_MAX_WINDOW
        windows = min_max_profile(codons, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        self.assertEqual(len(windows), 1)
        self.assertEqual(windows[0].percent, 100.0)
        self.assertEqual(windows[0].start_codon, 1)
        self.assertEqual(windows[0].end_codon, MIN_MAX_WINDOW)

    def test_all_lowest_preferred_codons_score_near_negative_100(self) -> None:
        codons = ["GCG"] * MIN_MAX_WINDOW
        windows = min_max_profile(codons, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        self.assertEqual(windows[0].percent, -100.0)

    def test_profile_is_empty_when_shorter_than_window(self) -> None:
        windows = min_max_profile(["GCT"] * (MIN_MAX_WINDOW - 1), PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        self.assertEqual(windows, [])

    def test_profile_slides_one_codon_at_a_time(self) -> None:
        codons = ["GCT"] * (MIN_MAX_WINDOW + 2)
        windows = min_max_profile(codons, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        self.assertEqual(len(windows), 3)
        self.assertEqual([window.start_codon for window in windows], [1, 2, 3])

    def test_window_with_no_comparable_codons_is_reported_not_dropped(self) -> None:
        codons = ["ATG"] * MIN_MAX_WINDOW  # Met has only one codon, no faster/slower choice
        windows = min_max_profile(codons, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        self.assertEqual(len(windows), 1)
        self.assertEqual(windows[0].start_codon, 1)
        self.assertEqual(windows[0].end_codon, MIN_MAX_WINDOW)
        self.assertIsNone(windows[0].percent)

    def test_start_codon_numbering_stays_contiguous_across_a_gap(self) -> None:
        codons = ["ATG"] * MIN_MAX_WINDOW + ["GCT"] * MIN_MAX_WINDOW
        windows = min_max_profile(codons, PUBLIC_PICHIA_PASTORIS_FRACTIONS)
        self.assertEqual([window.start_codon for window in windows], list(range(1, len(codons) - MIN_MAX_WINDOW + 2)))
        self.assertIsNone(windows[0].percent)
        self.assertIsNotNone(windows[-1].percent)


def _make_candidate(*, rank: int, cds: str, reference_cds: str) -> CdsCandidate:
    analysis = analyze_cds(cds)
    return CdsCandidate(
        rank=rank,
        generation_index=rank,
        source="test",
        cds=cds,
        codon_ids=[],
        analysis=analysis,
        quality=quality_summary(analysis),
        difference_from_reference=compare_cds(cds, reference_cds),
        codon_preference=codon_preference_stats(cds),
    )


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
                # NOTE: equal weights -> multinomial normalizes them, so this
                # draws UNIFORMLY and `temperature` has no effect here. The real
                # predictor uses softmax(logits / temperature), which concentrates
                # mass on preferred codons. This fake is therefore a worst case
                # for any gate that penalizes rare codons (e.g. the ADR-0001
                # avoidable-lowest check) -- do not use it to estimate how many
                # candidates the temperature-sampling strategy yields in practice.
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

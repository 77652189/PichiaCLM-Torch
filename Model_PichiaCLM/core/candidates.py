from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from itertools import combinations
import math
import random
import re
from typing import Iterable, Protocol

import torch

from .analysis import (
    AA_TO_CODONS,
    GLOBAL_GC_MAX,
    GLOBAL_GC_MIN,
    HOMOPOLYMER_MIN_LENGTH,
    PUBLIC_PICHIA_PASTORIS_FRACTIONS,
    SequenceAnalysisReport,
    analyze_cds,
    gc_percent,
)
from .biology import check_translation, normalize_motifs, split_codons
from .restriction import parse_custom_sites, scan_restriction_sites
from .schemas import PredictionResult
from .vocab import build_vocabularies


@dataclass(frozen=True)
class CandidateGenerationOptions:
    num_candidates: int = 10
    temperature: float = 0.8
    seed: int | None = None
    max_attempts: int | None = None
    strategy: str = "kazusa_diverse"
    pool_size: int | None = None
    min_difference_percent: float = 10.0
    max_difference_percent: float = 20.0
    subset_size: int | None = 5


@dataclass(frozen=True)
class CandidateQuality:
    status: str
    critical_issues: int
    warnings: int


@dataclass(frozen=True)
class CandidateDifference:
    bp_differences: int
    bp_difference_percent: float
    codon_differences: int
    codon_difference_percent: float


@dataclass(frozen=True)
class CandidatePairSimilarity:
    left_rank: int
    right_rank: int
    bp_similarity_percent: float
    bp_difference_percent: float
    codon_similarity_percent: float
    codon_difference_percent: float


@dataclass(frozen=True)
class CodonPreferenceStats:
    reference: str
    codon_count: int
    top_preferred_count: int
    second_preferred_count: int
    lowest_preferred_count: int
    avoidable_lowest_count: int
    top_preferred_percent: float
    second_preferred_percent: float
    lowest_preferred_percent: float
    avoidable_lowest_percent: float
    mean_fraction: float


@dataclass(frozen=True)
class PairwiseDiversity:
    comparisons: int
    min_bp_difference_percent: float | None
    mean_bp_difference_percent: float | None
    min_codon_difference_percent: float | None
    mean_codon_difference_percent: float | None


@dataclass(frozen=True)
class CandidateSubsetSelection:
    requested_size: int
    selected_size: int
    selected_ranks: list[int]
    method: str
    comparisons: int
    min_bp_similarity_percent: float | None
    mean_bp_similarity_percent: float | None
    max_bp_similarity_percent: float | None
    min_codon_similarity_percent: float | None
    mean_codon_similarity_percent: float | None
    max_codon_similarity_percent: float | None


@dataclass(frozen=True)
class CdsCandidate:
    rank: int
    generation_index: int
    source: str
    cds: str
    codon_ids: list[int]
    analysis: SequenceAnalysisReport
    quality: CandidateQuality
    difference_from_reference: CandidateDifference
    codon_preference: CodonPreferenceStats


@dataclass(frozen=True)
class CandidateSet:
    amino_acids: str
    reference_cds: str
    requested_candidates: int
    generated_candidates: int
    attempts: int
    exhausted: bool
    note: str | None
    pairwise_diversity: PairwiseDiversity
    pairwise_similarities: list[CandidatePairSimilarity]
    recommended_subset: CandidateSubsetSelection | None
    candidates: list[CdsCandidate]


@dataclass(frozen=True)
class CandidateDraft:
    generation_index: int
    cds: str
    codon_ids: list[int]
    changes: dict[int, str]
    difference_from_reference: CandidateDifference
    codon_preference: CodonPreferenceStats


@dataclass(frozen=True)
class LightweightRisk:
    gc_status: str
    homopolymers: int
    repeated_kmers: int
    motif_hits: int
    restriction_sites: int


class CandidatePredictor(Protocol):
    device: torch.device

    def predict(self, amino_acids: str, allow_unknown: bool = False) -> PredictionResult:
        ...

    def predict_sample(
        self,
        amino_acids: str,
        allow_unknown: bool = False,
        temperature: float = 0.8,
        generator: torch.Generator | None = None,
    ) -> PredictionResult:
        ...


def generate_cds_candidates(
    predictor: CandidatePredictor,
    amino_acids: str,
    *,
    options: CandidateGenerationOptions | None = None,
    allow_unknown: bool = False,
    motifs: Iterable[str] | None = None,
    custom_restriction_sites: Iterable[str] | None = None,
) -> CandidateSet:
    options = options or CandidateGenerationOptions()
    _validate_options(options)

    requested = options.num_candidates
    rng = random.Random(options.seed)
    reference = predictor.predict(amino_acids, allow_unknown=allow_unknown)
    if not _is_translation_consistent(reference):
        raise ValueError("Reference CDS does not translate back to the input amino acid sequence.")

    reference_analysis = analyze_cds(
        reference.cds,
        amino_acids=reference.amino_acids,
        motifs=motifs,
        custom_restriction_sites=custom_restriction_sites,
    )
    reference_quality = quality_summary(reference_analysis)
    reference_preference = codon_preference_stats(reference.cds)
    reference_risk = lightweight_risk(reference.cds, motifs, custom_restriction_sites)

    pool, attempts = _generate_candidate_pool(
        reference,
        options=options,
        rng=rng,
        reference_preference=reference_preference,
        reference_risk=reference_risk,
        motifs=motifs,
        custom_restriction_sites=custom_restriction_sites,
    )
    candidates = [
        _candidate_from_prediction(
            generation_index=1,
            source="reference",
            result=reference,
            reference_cds=reference.cds,
            analysis=reference_analysis,
        )
    ]

    selected_drafts: list[CandidateDraft] = []
    remaining_drafts = pool.copy()
    reference_codons = split_codons(reference.cds)
    while remaining_drafts and len(candidates) < requested:
        draft = _pop_most_diverse_draft(
            remaining_drafts,
            selected=selected_drafts,
            reference_codons=reference_codons,
        )
        candidate = _candidate_from_draft(
            draft,
            reference=reference,
            motifs=motifs,
            custom_restriction_sites=custom_restriction_sites,
        )
        if not _final_candidate_is_acceptable(candidate, reference_quality, reference_preference):
            continue
        selected_drafts.append(draft)
        candidates.append(candidate)

    candidates = [replace(candidate, rank=rank) for rank, candidate in enumerate(candidates, start=1)]
    pairwise_similarities = pairwise_similarity_rows(candidates)
    recommended_subset = select_low_similarity_subset(
        candidates,
        pairwise_similarities,
        subset_size=options.subset_size,
    )
    exhausted = len(candidates) < requested
    note = None
    if exhausted:
        note = (
            f"Requested {requested} unique candidates, but generated {len(candidates)} "
            f"after {attempts} lightweight attempts. The synonymous design space may be limited."
        )

    return CandidateSet(
        amino_acids=reference.amino_acids,
        reference_cds=reference.cds,
        requested_candidates=requested,
        generated_candidates=len(candidates),
        attempts=attempts,
        exhausted=exhausted,
        note=note,
        pairwise_diversity=pairwise_diversity_from_similarity(pairwise_similarities),
        pairwise_similarities=pairwise_similarities,
        recommended_subset=recommended_subset,
        candidates=candidates,
    )


def candidate_summary_rows(candidate_set: CandidateSet) -> list[dict[str, object]]:
    rows = []
    selected_ranks = set()
    if candidate_set.recommended_subset is not None:
        selected_ranks = set(candidate_set.recommended_subset.selected_ranks)
    for candidate in candidate_set.candidates:
        analysis = candidate.analysis
        similarity_values = [
            row.codon_similarity_percent
            for row in candidate_set.pairwise_similarities
            if row.left_rank == candidate.rank or row.right_rank == candidate.rank
        ]
        rows.append(
            {
                "rank": candidate.rank,
                "source": candidate.source,
                "recommended_subset": candidate.rank in selected_ranks,
                "aa_length": len(candidate_set.amino_acids),
                "cds_length": analysis.cds_length,
                "quality_status": candidate.quality.status,
                "critical_issues": candidate.quality.critical_issues,
                "warnings": candidate.quality.warnings,
                "gc_percent": analysis.gc_percent,
                "gc_status": analysis.gc_status,
                "cai_training": analysis.cai.training,
                "cai_public": analysis.cai.public,
                "kazusa_top_preferred_percent": candidate.codon_preference.top_preferred_percent,
                "kazusa_second_preferred_percent": candidate.codon_preference.second_preferred_percent,
                "kazusa_lowest_preferred_percent": candidate.codon_preference.lowest_preferred_percent,
                "kazusa_avoidable_lowest_percent": candidate.codon_preference.avoidable_lowest_percent,
                "kazusa_mean_fraction": candidate.codon_preference.mean_fraction,
                "bp_differences_from_reference": candidate.difference_from_reference.bp_differences,
                "bp_difference_percent": candidate.difference_from_reference.bp_difference_percent,
                "codon_differences_from_reference": candidate.difference_from_reference.codon_differences,
                "codon_difference_percent": candidate.difference_from_reference.codon_difference_percent,
                "mean_pairwise_codon_similarity_percent": _mean(similarity_values),
                "max_pairwise_codon_similarity_percent": max(similarity_values) if similarity_values else None,
                "restriction_sites": len(analysis.restriction_sites),
                "motif_hits": len(analysis.motif_hits),
                "local_gc_warnings": len(analysis.local_gc_outliers),
                "rare_codon_runs": len(analysis.rare_codon_runs),
                "homopolymers": len(analysis.homopolymers),
                "tandem_repeats": len(analysis.tandem_repeats),
                "repeated_kmers": len(analysis.repeated_kmers),
            }
        )
    return rows


def pairwise_similarity_rows(candidates: list[CdsCandidate]) -> list[CandidatePairSimilarity]:
    rows = []
    for left, right in combinations(candidates, 2):
        diff = compare_cds(left.cds, right.cds)
        rows.append(
            CandidatePairSimilarity(
                left_rank=left.rank,
                right_rank=right.rank,
                bp_similarity_percent=round(100.0 - diff.bp_difference_percent, 2),
                bp_difference_percent=diff.bp_difference_percent,
                codon_similarity_percent=round(100.0 - diff.codon_difference_percent, 2),
                codon_difference_percent=diff.codon_difference_percent,
            )
        )
    return rows


def pairwise_diversity(candidates: list[CdsCandidate]) -> PairwiseDiversity:
    return pairwise_diversity_from_similarity(pairwise_similarity_rows(candidates))


def pairwise_diversity_from_similarity(rows: list[CandidatePairSimilarity]) -> PairwiseDiversity:
    if not rows:
        return PairwiseDiversity(
            comparisons=0,
            min_bp_difference_percent=None,
            mean_bp_difference_percent=None,
            min_codon_difference_percent=None,
            mean_codon_difference_percent=None,
        )
    bp_values = [row.bp_difference_percent for row in rows]
    codon_values = [row.codon_difference_percent for row in rows]
    return PairwiseDiversity(
        comparisons=len(rows),
        min_bp_difference_percent=round(min(bp_values), 2),
        mean_bp_difference_percent=round(sum(bp_values) / len(bp_values), 2),
        min_codon_difference_percent=round(min(codon_values), 2),
        mean_codon_difference_percent=round(sum(codon_values) / len(codon_values), 2),
    )


def select_low_similarity_subset(
    candidates: list[CdsCandidate],
    pairwise_rows: list[CandidatePairSimilarity],
    *,
    subset_size: int | None = 5,
    max_exact_combinations: int = 20_000,
) -> CandidateSubsetSelection | None:
    if not candidates or subset_size is None:
        return None
    requested_size = max(1, subset_size)
    selected_size = min(requested_size, len(candidates))
    ranks = [candidate.rank for candidate in candidates]
    pair_lookup = _pairwise_lookup(pairwise_rows)

    if selected_size == len(ranks):
        return _subset_selection_from_ranks(
            ranks,
            requested_size=requested_size,
            method="all_candidates",
            pair_lookup=pair_lookup,
        )

    combination_count = math.comb(len(ranks), selected_size)
    if combination_count <= max_exact_combinations:
        best_ranks = min(
            combinations(ranks, selected_size),
            key=lambda item: _subset_score(item, pair_lookup),
        )
        method = "exact"
    else:
        best_ranks = _greedy_low_similarity_subset(ranks, selected_size, pair_lookup)
        method = "greedy"

    return _subset_selection_from_ranks(
        list(best_ranks),
        requested_size=requested_size,
        method=method,
        pair_lookup=pair_lookup,
    )


def _pairwise_lookup(
    rows: list[CandidatePairSimilarity],
) -> dict[tuple[int, int], CandidatePairSimilarity]:
    return {
        tuple(sorted((row.left_rank, row.right_rank))): row
        for row in rows
    }


def _subset_score(
    ranks: Iterable[int],
    pair_lookup: dict[tuple[int, int], CandidatePairSimilarity],
) -> tuple[float, float, float, int]:
    pair_rows = _subset_pair_rows(ranks, pair_lookup)
    if not pair_rows:
        return (0.0, 0.0, 0.0, sum(ranks))
    codon_similarities = [row.codon_similarity_percent for row in pair_rows]
    return (
        max(codon_similarities),
        sum(codon_similarities) / len(codon_similarities),
        -min(row.codon_difference_percent for row in pair_rows),
        sum(ranks),
    )


def _greedy_low_similarity_subset(
    ranks: list[int],
    selected_size: int,
    pair_lookup: dict[tuple[int, int], CandidatePairSimilarity],
) -> list[int]:
    if selected_size <= 1:
        return ranks[:1]
    selected = list(
        min(
            combinations(ranks, 2),
            key=lambda item: _subset_score(item, pair_lookup),
        )
    )
    while len(selected) < selected_size:
        remaining = [rank for rank in ranks if rank not in selected]
        next_rank = min(
            remaining,
            key=lambda rank: _subset_score([*selected, rank], pair_lookup),
        )
        selected.append(next_rank)
    return selected


def _subset_selection_from_ranks(
    ranks: list[int] | tuple[int, ...],
    *,
    requested_size: int,
    method: str,
    pair_lookup: dict[tuple[int, int], CandidatePairSimilarity],
) -> CandidateSubsetSelection:
    selected_ranks = sorted(ranks)
    pair_rows = _subset_pair_rows(selected_ranks, pair_lookup)
    bp_similarities = [row.bp_similarity_percent for row in pair_rows]
    codon_similarities = [row.codon_similarity_percent for row in pair_rows]
    return CandidateSubsetSelection(
        requested_size=requested_size,
        selected_size=len(selected_ranks),
        selected_ranks=selected_ranks,
        method=method,
        comparisons=len(pair_rows),
        min_bp_similarity_percent=min(bp_similarities) if bp_similarities else None,
        mean_bp_similarity_percent=_mean(bp_similarities),
        max_bp_similarity_percent=max(bp_similarities) if bp_similarities else None,
        min_codon_similarity_percent=min(codon_similarities) if codon_similarities else None,
        mean_codon_similarity_percent=_mean(codon_similarities),
        max_codon_similarity_percent=max(codon_similarities) if codon_similarities else None,
    )


def _subset_pair_rows(
    ranks: Iterable[int],
    pair_lookup: dict[tuple[int, int], CandidatePairSimilarity],
) -> list[CandidatePairSimilarity]:
    rows = []
    for left, right in combinations(sorted(ranks), 2):
        pair = tuple(sorted((left, right)))
        if pair in pair_lookup:
            rows.append(pair_lookup[pair])
    return rows


def compare_cds(cds: str, reference_cds: str) -> CandidateDifference:
    bp_differences = _hamming_with_length(cds, reference_cds)
    bp_total = max(len(cds), len(reference_cds), 1)
    codons = split_codons(cds)
    reference_codons = split_codons(reference_cds)
    codon_differences = _hamming_with_length(codons, reference_codons)
    codon_total = max(len(codons), len(reference_codons), 1)
    return CandidateDifference(
        bp_differences=bp_differences,
        bp_difference_percent=round(bp_differences / bp_total * 100, 2),
        codon_differences=codon_differences,
        codon_difference_percent=round(codon_differences / codon_total * 100, 2),
    )


def quality_summary(analysis: SequenceAnalysisReport) -> CandidateQuality:
    critical = 0
    if not analysis.valid_dna:
        critical += 1
    if not analysis.length_multiple_of_three:
        critical += 1
    if analysis.translation_matches_input is False:
        critical += 1
    critical += len(analysis.internal_stop_codons)

    warnings = len(analysis.sequence_warnings)
    warnings += 1 if analysis.gc_status != "ok" else 0
    warnings += len(analysis.local_gc_outliers)
    warnings += len(analysis.restriction_sites)
    warnings += len(analysis.motif_hits)
    warnings += len(analysis.rare_codon_runs)
    warnings += len(analysis.homopolymers)
    warnings += len(analysis.tandem_repeats)
    warnings += len(analysis.repeated_kmers)

    if critical:
        status = "fail"
    elif warnings:
        status = "warning"
    else:
        status = "pass"
    return CandidateQuality(status=status, critical_issues=critical, warnings=warnings)


def codon_preference_stats(cds: str) -> CodonPreferenceStats:
    codons = split_codons(cds)
    top_count = 0
    second_count = 0
    lowest_count = 0
    avoidable_lowest_count = 0
    fractions = []

    for codon in codons:
        fraction = PUBLIC_PICHIA_PASTORIS_FRACTIONS.get(codon, 0.0)
        fractions.append(fraction)
        ranked = _ranked_synonymous_codons(codon)
        if not ranked:
            continue
        rank = ranked.index(codon) if codon in ranked else -1
        if rank == 0:
            top_count += 1
        if rank == 1:
            second_count += 1
        if rank == len(ranked) - 1 and len(ranked) > 1:
            lowest_count += 1
            avoidable_lowest_count += 1

    total = max(len(codons), 1)
    return CodonPreferenceStats(
        reference="Kazusa Pichia pastoris taxon 4922",
        codon_count=len(codons),
        top_preferred_count=top_count,
        second_preferred_count=second_count,
        lowest_preferred_count=lowest_count,
        avoidable_lowest_count=avoidable_lowest_count,
        top_preferred_percent=round(top_count / total * 100, 2),
        second_preferred_percent=round(second_count / total * 100, 2),
        lowest_preferred_percent=round(lowest_count / total * 100, 2),
        avoidable_lowest_percent=round(avoidable_lowest_count / total * 100, 2),
        mean_fraction=round(sum(fractions) / total, 4),
    )


def lightweight_risk(
    cds: str,
    motifs: Iterable[str] | None = None,
    custom_restriction_sites: Iterable[str] | None = None,
) -> LightweightRisk:
    global_gc = gc_percent(cds)
    if global_gc < GLOBAL_GC_MIN:
        gc_status = "low"
    elif global_gc > GLOBAL_GC_MAX:
        gc_status = "high"
    else:
        gc_status = "ok"
    return LightweightRisk(
        gc_status=gc_status,
        homopolymers=_count_homopolymers(cds),
        repeated_kmers=_count_repeated_kmers(cds),
        motif_hits=_count_motif_hits(cds, motifs),
        restriction_sites=len(
            scan_restriction_sites(
                cds,
                include_defaults=True,
                custom_sites=parse_custom_sites(custom_restriction_sites),
            )
        ),
    )


def _generate_candidate_pool(
    reference: PredictionResult,
    *,
    options: CandidateGenerationOptions,
    rng: random.Random,
    reference_preference: CodonPreferenceStats,
    reference_risk: LightweightRisk,
    motifs: Iterable[str] | None,
    custom_restriction_sites: Iterable[str] | None,
) -> tuple[list[CandidateDraft], int]:
    reference_codons = split_codons(reference.cds)
    replacements_by_position = {
        index: _replacement_codons(codon)
        for index, codon in enumerate(reference_codons)
        if _replacement_codons(codon)
    }
    mutable_positions = list(replacements_by_position)
    min_changes, max_changes = _change_count_bounds(len(reference_codons), options)
    if not mutable_positions or max_changes < 1:
        return [], 1

    min_changes = min(min_changes, len(mutable_positions))
    max_changes = min(max_changes, len(mutable_positions))
    if min_changes > max_changes:
        return [], 1

    pool_size = options.pool_size or min(120, max(20, (options.num_candidates - 1) * 6))
    max_attempts = options.max_attempts or max(pool_size * 4, 60)
    stale_limit = max(pool_size * 2, 40)
    pool: list[CandidateDraft] = []
    seen_signatures: set[tuple[tuple[int, str], ...]] = set()
    attempts = 1
    stale_attempts = 0

    while len(pool) < pool_size and attempts < max_attempts and stale_attempts < stale_limit:
        attempts += 1
        change_count = rng.randint(min_changes, max_changes)
        positions = rng.sample(mutable_positions, change_count)
        changes = {}
        for position in positions:
            alternatives = replacements_by_position[position]
            changes[position] = _weighted_choice(alternatives, rng, temperature=options.temperature)
        signature = tuple(sorted(changes.items()))
        if signature in seen_signatures:
            stale_attempts += 1
            continue
        seen_signatures.add(signature)

        draft = _draft_from_changes(reference, reference_codons, changes, attempts)
        if not _draft_remains_lightweight_acceptable(
            draft,
            reference=reference,
            reference_preference=reference_preference,
            reference_risk=reference_risk,
            options=options,
            motifs=motifs,
            custom_restriction_sites=custom_restriction_sites,
        ):
            stale_attempts += 1
            continue
        pool.append(draft)
        stale_attempts = 0

    return pool, attempts


def _draft_from_changes(
    reference: PredictionResult,
    reference_codons: list[str],
    changes: dict[int, str],
    generation_index: int,
) -> CandidateDraft:
    candidate_codons = reference_codons.copy()
    for position, codon in changes.items():
        candidate_codons[position] = codon
    cds = "".join(candidate_codons)
    return CandidateDraft(
        generation_index=generation_index,
        cds=cds,
        codon_ids=_codon_ids_from_reference(reference.codon_ids, reference_codons, candidate_codons),
        changes=changes,
        difference_from_reference=compare_cds(cds, reference.cds),
        codon_preference=codon_preference_stats(cds),
    )


def _draft_remains_lightweight_acceptable(
    draft: CandidateDraft,
    *,
    reference: PredictionResult,
    reference_preference: CodonPreferenceStats,
    reference_risk: LightweightRisk,
    options: CandidateGenerationOptions,
    motifs: Iterable[str] | None,
    custom_restriction_sites: Iterable[str] | None,
) -> bool:
    if not _is_translation_consistent(
        PredictionResult(
            amino_acids=reference.amino_acids,
            cds=draft.cds,
            codon_ids=draft.codon_ids,
            device=reference.device,
        )
    ):
        return False
    diff = draft.difference_from_reference.codon_difference_percent
    if diff < options.min_difference_percent or diff > options.max_difference_percent:
        return False
    if draft.codon_preference.avoidable_lowest_count > reference_preference.avoidable_lowest_count:
        return False

    risk = lightweight_risk(draft.cds, motifs, custom_restriction_sites)
    if _gc_risk_score(risk.gc_status) > _gc_risk_score(reference_risk.gc_status):
        return False
    if risk.homopolymers > reference_risk.homopolymers:
        return False
    if risk.repeated_kmers > reference_risk.repeated_kmers:
        return False
    if risk.motif_hits > reference_risk.motif_hits:
        return False
    if risk.restriction_sites > reference_risk.restriction_sites:
        return False
    return True


def _select_diverse_drafts(
    pool: list[CandidateDraft],
    *,
    target: int,
    reference_cds: str,
) -> list[CandidateDraft]:
    if target <= 0:
        return []
    reference_codons = split_codons(reference_cds)
    remaining = pool.copy()
    selected: list[CandidateDraft] = []
    while remaining and len(selected) < target:
        best = max(
            remaining,
            key=lambda draft: _diversity_score(draft, selected, reference_codons),
        )
        selected.append(best)
        remaining.remove(best)
    return selected


def _pop_most_diverse_draft(
    remaining: list[CandidateDraft],
    *,
    selected: list[CandidateDraft],
    reference_codons: list[str],
) -> CandidateDraft:
    best = max(
        remaining,
        key=lambda draft: _diversity_score(draft, selected, reference_codons),
    )
    remaining.remove(best)
    return best


def _diversity_score(
    draft: CandidateDraft,
    selected: list[CandidateDraft],
    reference_codons: list[str],
) -> tuple[float, float, float, float, int]:
    if selected:
        min_diversity = min(_draft_pair_difference_percent(draft, item, reference_codons) for item in selected)
    else:
        min_diversity = draft.difference_from_reference.codon_difference_percent
    return (
        min_diversity,
        draft.difference_from_reference.codon_difference_percent,
        draft.codon_preference.second_preferred_percent,
        draft.codon_preference.mean_fraction,
        -draft.generation_index,
    )


def _draft_pair_difference_percent(
    left: CandidateDraft,
    right: CandidateDraft,
    reference_codons: list[str],
) -> float:
    positions = set(left.changes) | set(right.changes)
    if not reference_codons:
        return 0.0
    differences = 0
    for position in positions:
        left_codon = left.changes.get(position, reference_codons[position])
        right_codon = right.changes.get(position, reference_codons[position])
        if left_codon != right_codon:
            differences += 1
    return round(differences / len(reference_codons) * 100, 2)


def _candidate_from_draft(
    draft: CandidateDraft,
    *,
    reference: PredictionResult,
    motifs: Iterable[str] | None,
    custom_restriction_sites: Iterable[str] | None,
) -> CdsCandidate:
    return _candidate_from_prediction(
        generation_index=draft.generation_index,
        source="kazusa_diverse",
        result=PredictionResult(
            amino_acids=reference.amino_acids,
            cds=draft.cds,
            codon_ids=draft.codon_ids,
            device=reference.device,
        ),
        reference_cds=reference.cds,
        analysis=None,
        motifs=motifs,
        custom_restriction_sites=custom_restriction_sites,
    )


def _candidate_from_prediction(
    *,
    generation_index: int,
    source: str,
    result: PredictionResult,
    reference_cds: str,
    analysis: SequenceAnalysisReport | None = None,
    motifs: Iterable[str] | None = None,
    custom_restriction_sites: Iterable[str] | None = None,
) -> CdsCandidate:
    analysis = analysis or analyze_cds(
        result.cds,
        amino_acids=result.amino_acids,
        motifs=motifs,
        custom_restriction_sites=custom_restriction_sites,
    )
    return CdsCandidate(
        rank=0,
        generation_index=generation_index,
        source=source,
        cds=result.cds,
        codon_ids=result.codon_ids,
        analysis=analysis,
        quality=quality_summary(analysis),
        difference_from_reference=compare_cds(result.cds, reference_cds),
        codon_preference=codon_preference_stats(result.cds),
    )


def _final_candidate_is_acceptable(
    candidate: CdsCandidate,
    reference_quality: CandidateQuality,
    reference_preference: CodonPreferenceStats,
) -> bool:
    if candidate.quality.critical_issues:
        return False
    if candidate.quality.warnings > reference_quality.warnings:
        return False
    if candidate.codon_preference.avoidable_lowest_count > reference_preference.avoidable_lowest_count:
        return False
    return True


def _is_translation_consistent(result: PredictionResult) -> bool:
    analysis = check_translation(result.cds, amino_acids=result.amino_acids)
    return (
        analysis.valid_dna
        and analysis.length_multiple_of_three
        and analysis.translation_matches_input is True
        and not analysis.internal_stop_codons
    )


def _replacement_codons(codon: str) -> list[str]:
    ranked = _ranked_synonymous_codons(codon)
    if len(ranked) <= 1 or codon not in ranked:
        return []
    lowest = ranked[-1]
    preferred = []
    if len(ranked) > 1 and ranked[1] != codon and ranked[1] != lowest:
        preferred.append(ranked[1])
    for item in ranked:
        if item != codon and item != lowest and item not in preferred:
            preferred.append(item)
    return preferred


def _ranked_synonymous_codons(codon: str) -> list[str]:
    aa = None
    for candidate_aa, codons in AA_TO_CODONS.items():
        if codon in codons:
            aa = candidate_aa
            break
    if aa is None:
        return []
    return sorted(
        AA_TO_CODONS.get(aa, []),
        key=lambda item: (PUBLIC_PICHIA_PASTORIS_FRACTIONS.get(item, 0.0), item),
        reverse=True,
    )


def _codon_ids_from_reference(
    reference_ids: list[int],
    reference_codons: list[str],
    candidate_codons: list[str],
) -> list[int]:
    _, codon_by_id, _ = build_vocabularies()
    id_by_codon = {codon: codon_id for codon_id, codon in codon_by_id.items()}
    return [
        reference_id if reference_codon == candidate_codon else id_by_codon.get(candidate_codon, 0)
        for reference_id, reference_codon, candidate_codon in zip(reference_ids, reference_codons, candidate_codons)
    ]


def _change_count_bounds(codon_count: int, options: CandidateGenerationOptions) -> tuple[int, int]:
    min_changes = math.ceil(codon_count * options.min_difference_percent / 100)
    max_changes = math.floor(codon_count * options.max_difference_percent / 100)
    return max(1, min_changes), max_changes


def _weighted_choice(items: list[str], rng: random.Random, *, temperature: float) -> str:
    base_weights = [len(items) - index for index in range(len(items))]
    weights = [weight ** (1.0 / temperature) for weight in base_weights]
    return rng.choices(items, weights=weights, k=1)[0]


def _gc_risk_score(status: str) -> int:
    return 0 if status == "ok" else 1


def _count_homopolymers(cds: str) -> int:
    return len(re.findall(r"([ATGC])\1{%d,}" % (HOMOPOLYMER_MIN_LENGTH - 1), cds))


def _count_repeated_kmers(cds: str, kmer_size: int = 12, min_count: int = 3) -> int:
    if len(cds) < kmer_size:
        return 0
    counts = Counter(cds[index : index + kmer_size] for index in range(0, len(cds) - kmer_size + 1))
    return sum(1 for count in counts.values() if count >= min_count)


def _count_motif_hits(cds: str, motifs: Iterable[str] | None) -> int:
    hits = 0
    for motif in normalize_motifs(motifs):
        start = 0
        while True:
            index = cds.find(motif, start)
            if index == -1:
                break
            hits += 1
            start = index + 1
    return hits


def _validate_options(options: CandidateGenerationOptions) -> None:
    if options.num_candidates < 1:
        raise ValueError("num_candidates must be at least 1.")
    if options.temperature <= 0:
        raise ValueError("temperature must be greater than 0.")
    if options.max_attempts is not None and options.max_attempts < options.num_candidates:
        raise ValueError("max_attempts must be greater than or equal to num_candidates.")
    if options.pool_size is not None and options.pool_size < options.num_candidates:
        raise ValueError("pool_size must be greater than or equal to num_candidates.")
    if options.subset_size is not None and options.subset_size < 1:
        raise ValueError("subset_size must be at least 1.")
    if options.min_difference_percent < 0:
        raise ValueError("min_difference_percent must be non-negative.")
    if options.max_difference_percent <= 0:
        raise ValueError("max_difference_percent must be greater than 0.")
    if options.min_difference_percent > options.max_difference_percent:
        raise ValueError("min_difference_percent must not exceed max_difference_percent.")


def _hamming_with_length(left, right) -> int:
    shared = sum(1 for left_item, right_item in zip(left, right) if left_item != right_item)
    return shared + abs(len(left) - len(right))


def _mean(values: list[float]) -> float | None:
    if not values:
        return None
    return round(sum(values) / len(values), 2)

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations
import random
from typing import Iterable, Protocol

import torch

from .analysis import AA_TO_CODONS, PUBLIC_PICHIA_PASTORIS_FRACTIONS, SequenceAnalysisReport, analyze_cds
from .biology import check_translation, split_codons
from .schemas import PredictionResult
from .vocab import build_vocabularies


@dataclass(frozen=True)
class CandidateGenerationOptions:
    num_candidates: int = 10
    temperature: float = 0.8
    seed: int | None = None
    max_attempts: int | None = None
    strategy: str = "kazusa_constrained"


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
    candidates: list[CdsCandidate]


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
    max_attempts = options.max_attempts or max(requested * 12, 60)
    rng = random.Random(options.seed)

    raw_results: list[tuple[int, str, PredictionResult]] = []
    seen_cds: set[str] = set()
    attempts = 1

    reference = predictor.predict(amino_acids, allow_unknown=allow_unknown)
    if _is_translation_consistent(reference):
        raw_results.append((1, "reference", reference))
        seen_cds.add(reference.cds)

    reference_analysis = analyze_cds(
        reference.cds,
        amino_acids=reference.amino_acids,
        motifs=motifs,
        custom_restriction_sites=custom_restriction_sites,
    )
    reference_quality = quality_summary(reference_analysis)
    reference_preference = codon_preference_stats(reference.cds)

    while len(raw_results) < requested and attempts < max_attempts:
        attempts += 1
        candidate = _build_kazusa_constrained_variant(
            reference,
            attempt=attempts,
            rng=rng,
            reference_quality=reference_quality,
            max_avoidable_lowest_count=reference_preference.avoidable_lowest_count,
            motifs=motifs,
            custom_restriction_sites=custom_restriction_sites,
        )
        if candidate is None or candidate.cds in seen_cds:
            continue
        if not _is_translation_consistent(candidate):
            continue
        raw_results.append((attempts, "kazusa_constrained", candidate))
        seen_cds.add(candidate.cds)

    candidates = _build_candidates(
        raw_results,
        reference_cds=reference.cds,
        motifs=motifs,
        custom_restriction_sites=custom_restriction_sites,
    )
    exhausted = len(candidates) < requested
    note = None
    if exhausted:
        note = (
            f"Requested {requested} unique candidates, but generated {len(candidates)} "
            f"after {attempts} attempts. The synonymous design space may be limited."
        )

    return CandidateSet(
        amino_acids=reference.amino_acids,
        reference_cds=reference.cds,
        requested_candidates=requested,
        generated_candidates=len(candidates),
        attempts=attempts,
        exhausted=exhausted,
        note=note,
        pairwise_diversity=pairwise_diversity(candidates),
        candidates=candidates,
    )


def candidate_summary_rows(candidate_set: CandidateSet) -> list[dict[str, object]]:
    rows = []
    for candidate in candidate_set.candidates:
        analysis = candidate.analysis
        rows.append(
            {
                "rank": candidate.rank,
                "source": candidate.source,
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


def pairwise_diversity(candidates: list[CdsCandidate]) -> PairwiseDiversity:
    bp_values = []
    codon_values = []
    for left, right in combinations(candidates, 2):
        diff = compare_cds(left.cds, right.cds)
        bp_values.append(diff.bp_difference_percent)
        codon_values.append(diff.codon_difference_percent)

    if not bp_values:
        return PairwiseDiversity(
            comparisons=0,
            min_bp_difference_percent=None,
            mean_bp_difference_percent=None,
            min_codon_difference_percent=None,
            mean_codon_difference_percent=None,
        )

    return PairwiseDiversity(
        comparisons=len(bp_values),
        min_bp_difference_percent=round(min(bp_values), 2),
        mean_bp_difference_percent=round(sum(bp_values) / len(bp_values), 2),
        min_codon_difference_percent=round(min(codon_values), 2),
        mean_codon_difference_percent=round(sum(codon_values) / len(codon_values), 2),
    )


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


def _build_candidates(
    raw_results: list[tuple[int, str, PredictionResult]],
    *,
    reference_cds: str,
    motifs: Iterable[str] | None,
    custom_restriction_sites: Iterable[str] | None,
) -> list[CdsCandidate]:
    candidates = []
    for generation_index, source, result in raw_results:
        analysis = analyze_cds(
            result.cds,
            amino_acids=result.amino_acids,
            motifs=motifs,
            custom_restriction_sites=custom_restriction_sites,
        )
        candidates.append(
            CdsCandidate(
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
        )

    reference_candidates = [candidate for candidate in candidates if candidate.source == "reference"]
    sampled_candidates = [candidate for candidate in candidates if candidate.source != "reference"]
    sampled_candidates.sort(
        key=lambda item: (
            item.quality.critical_issues,
            item.quality.warnings,
            -item.difference_from_reference.codon_difference_percent,
            item.generation_index,
        )
    )
    candidates = reference_candidates + sampled_candidates
    return [replace(candidate, rank=rank) for rank, candidate in enumerate(candidates, start=1)]


def _is_translation_consistent(result: PredictionResult) -> bool:
    analysis = check_translation(result.cds, amino_acids=result.amino_acids)
    return (
        analysis.valid_dna
        and analysis.length_multiple_of_three
        and analysis.translation_matches_input is True
        and not analysis.internal_stop_codons
    )


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


def _build_kazusa_constrained_variant(
    reference: PredictionResult,
    *,
    attempt: int,
    rng: random.Random,
    reference_quality: CandidateQuality,
    max_avoidable_lowest_count: int,
    motifs: Iterable[str] | None,
    custom_restriction_sites: Iterable[str] | None,
) -> PredictionResult | None:
    codons = split_codons(reference.cds)
    mutable_positions = [
        index
        for index, codon in enumerate(codons)
        if _replacement_codons(codon)
    ]
    if not mutable_positions:
        return None

    rng.shuffle(mutable_positions)
    replacement_fraction = 0.10 + (attempt % 4) * 0.03
    target_replacements = max(1, min(len(mutable_positions), round(len(mutable_positions) * replacement_fraction)))
    trial_codons = codons.copy()
    replacements = 0

    for position in mutable_positions:
        old_codon = trial_codons[position]
        alternatives = _replacement_codons(old_codon)
        if not alternatives:
            continue
        rotation = attempt % len(alternatives)
        alternatives = alternatives[rotation:] + alternatives[:rotation]
        for new_codon in alternatives:
            next_codons = trial_codons.copy()
            next_codons[position] = new_codon
            trial_cds = "".join(next_codons)
            if compare_cds(trial_cds, reference.cds).codon_difference_percent > 20.0:
                continue
            if not _candidate_remains_acceptable(
                trial_cds,
                reference.amino_acids,
                reference_quality,
                max_avoidable_lowest_count,
                motifs,
                custom_restriction_sites,
            ):
                continue
            trial_codons = next_codons
            replacements += 1
            break
        if replacements >= target_replacements:
            break

    if replacements == 0:
        return None
    candidate_cds = "".join(trial_codons)
    return PredictionResult(
        amino_acids=reference.amino_acids,
        cds=candidate_cds,
        codon_ids=_codon_ids_from_reference(reference.codon_ids, codons, trial_codons),
        device=reference.device,
    )


def _candidate_remains_acceptable(
    cds: str,
    amino_acids: str,
    reference_quality: CandidateQuality,
    max_avoidable_lowest_count: int,
    motifs: Iterable[str] | None,
    custom_restriction_sites: Iterable[str] | None,
) -> bool:
    translation = check_translation(cds, amino_acids=amino_acids)
    if (
        not translation.valid_dna
        or not translation.length_multiple_of_three
        or translation.translation_matches_input is not True
        or translation.internal_stop_codons
    ):
        return False
    preference = codon_preference_stats(cds)
    if preference.avoidable_lowest_count > max_avoidable_lowest_count:
        return False
    analysis = analyze_cds(
        cds,
        amino_acids=amino_acids,
        motifs=motifs,
        custom_restriction_sites=custom_restriction_sites,
    )
    quality = quality_summary(analysis)
    if quality.critical_issues:
        return False
    return quality.warnings <= reference_quality.warnings


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


def _codon_ids_from_reference(reference_ids: list[int], reference_codons: list[str], candidate_codons: list[str]) -> list[int]:
    _, codon_by_id, _ = build_vocabularies()
    id_by_codon = {codon: codon_id for codon_id, codon in codon_by_id.items()}
    return [
        reference_id if reference_codon == candidate_codon else id_by_codon.get(candidate_codon, 0)
        for reference_id, reference_codon, candidate_codon in zip(reference_ids, reference_codons, candidate_codons)
    ]


def _validate_options(options: CandidateGenerationOptions) -> None:
    if options.num_candidates < 1:
        raise ValueError("num_candidates must be at least 1.")
    if options.temperature <= 0:
        raise ValueError("temperature must be greater than 0.")
    if options.max_attempts is not None and options.max_attempts < options.num_candidates:
        raise ValueError("max_attempts must be greater than or equal to num_candidates.")


def _hamming_with_length(left, right) -> int:
    shared = sum(1 for left_item, right_item in zip(left, right) if left_item != right_item)
    return shared + abs(len(left) - len(right))

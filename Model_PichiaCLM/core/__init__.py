"""Core PichiaCLM inference package.

This package intentionally has no dependency on FastAPI, Streamlit, or CLI code.
"""

from .analysis import SequenceAnalysisReport, analyze_cds
from .biology import TranslationCheck, check_translation, translate_cds
from .candidates import (
    CandidateGenerationOptions,
    CandidateSet,
    CandidateSubsetSelection,
    CdsCandidate,
    generate_cds_candidates,
    select_low_similarity_subset,
)
from .fasta import FastaRecord, format_fasta, parse_fasta
from .fusion import FusionComparison, compare_signal_fusion
from .postprocess import PostprocessResult, conservative_postprocess
from .predictor import PichiaCLMPredictor, batch_predict
from .restriction import RestrictionSite, scan_restriction_sites
from .schemas import PredictionResult

__all__ = [
    "FastaRecord",
    "CandidateGenerationOptions",
    "CandidateSet",
    "CandidateSubsetSelection",
    "CdsCandidate",
    "FusionComparison",
    "PichiaCLMPredictor",
    "PostprocessResult",
    "PredictionResult",
    "RestrictionSite",
    "SequenceAnalysisReport",
    "TranslationCheck",
    "analyze_cds",
    "batch_predict",
    "check_translation",
    "compare_signal_fusion",
    "conservative_postprocess",
    "format_fasta",
    "generate_cds_candidates",
    "parse_fasta",
    "scan_restriction_sites",
    "select_low_similarity_subset",
    "translate_cds",
]

"""Backward-compatible imports for the old inference module path."""

from .core.config import DEFAULT_WEIGHTS_PATH
from .core.model import MultiTaskSeq2Seq
from .core.predictor import PichiaCLMPredictor, batch_predict
from .core.schemas import PredictionResult
from .core.vocab import (
    AA_EOS_IDX,
    AA_UNK_IDX,
    CDS_SOS_IDX,
    PAD_IDX,
    build_vocabularies,
    normalize_amino_acids,
)

__all__ = [
    "AA_EOS_IDX",
    "AA_UNK_IDX",
    "CDS_SOS_IDX",
    "DEFAULT_WEIGHTS_PATH",
    "MultiTaskSeq2Seq",
    "PAD_IDX",
    "PichiaCLMPredictor",
    "PredictionResult",
    "batch_predict",
    "build_vocabularies",
    "normalize_amino_acids",
]

"""Backward-compatible FastAPI module.

Prefer: uvicorn Model_PichiaCLM.interfaces.api:app
"""

from .interfaces.api import PredictRequest, PredictResponse, app, get_predictor, health, predict

__all__ = ["PredictRequest", "PredictResponse", "app", "get_predictor", "health", "predict"]

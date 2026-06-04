from __future__ import annotations

import os
from dataclasses import asdict
from functools import lru_cache

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from Model_PichiaCLM.core.config import DEFAULT_WEIGHTS_PATH
from Model_PichiaCLM.core.predictor import PichiaCLMPredictor


class PredictRequest(BaseModel):
    amino_acids: str = Field(..., examples=["MSTNPKPQR"])
    allow_unknown: bool = False


class PredictResponse(BaseModel):
    amino_acids: str
    cds: str
    codon_ids: list[int]
    device: str


@lru_cache(maxsize=1)
def get_predictor() -> PichiaCLMPredictor:
    weights_path = os.environ.get("PICHIA_CLM_WEIGHTS", str(DEFAULT_WEIGHTS_PATH))
    device = os.environ.get("PICHIA_CLM_DEVICE") or None
    return PichiaCLMPredictor(weights_path=weights_path, device=device)


app = FastAPI(title="PichiaCLM Torch API", version="0.1.0")


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/predict", response_model=PredictResponse)
def predict(request: PredictRequest) -> dict[str, object]:
    try:
        result = get_predictor().predict(
            request.amino_acids,
            allow_unknown=request.allow_unknown,
        )
    except (FileNotFoundError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return asdict(result)

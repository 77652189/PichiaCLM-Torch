# PichiaCLM Torch Deployment

This repository ships research notebooks plus trained PyTorch weights. The deployment code is split so the model logic stays independent from each user interface.

## Architecture

```text
Model_PichiaCLM/
  core/
    config.py       # default model paths
    model.py        # PyTorch model architecture
    predictor.py    # model loading and AA-to-CDS prediction
    schemas.py      # shared result dataclasses
    vocab.py        # AA/CDS vocabularies and biological masks
  interfaces/
    cli.py          # command-line interface
    api.py          # FastAPI interface
    streamlit_app.py # Streamlit interface
```

`core` does not import FastAPI, Streamlit, or CLI code. CLI, API, and Streamlit are thin adapters over the same `PichiaCLMPredictor`.

Legacy module paths still work:

```powershell
python -m Model_PichiaCLM.cli --aa MSTNPKPQR
uvicorn Model_PichiaCLM.api:app --host 127.0.0.1 --port 8000
```

Prefer the new explicit paths:

```powershell
python -m Model_PichiaCLM.interfaces.cli --aa MSTNPKPQR
uvicorn Model_PichiaCLM.interfaces.api:app --host 127.0.0.1 --port 8000
streamlit run Model_PichiaCLM/interfaces/streamlit_app.py
```

## Install

Core inference only:

```powershell
pip install -r requirements-core.txt
```

FastAPI service:

```powershell
pip install -r requirements-api.txt
```

Streamlit UI:

```powershell
pip install -r requirements-streamlit.txt
```

Everything:

```powershell
pip install -r requirements-deploy.txt
```

## CLI

```powershell
python -m Model_PichiaCLM.interfaces.cli --aa MSTNPKPQR
python -m Model_PichiaCLM.interfaces.cli --aa MSTNPKPQR --json
```

## API

```powershell
uvicorn Model_PichiaCLM.interfaces.api:app --host 127.0.0.1 --port 8000
```

Call the API:

```powershell
Invoke-RestMethod -Method Post `
  -Uri http://127.0.0.1:8000/predict `
  -ContentType application/json `
  -Body '{"amino_acids":"MSTNPKPQR"}'
```

## Streamlit

Direct model mode:

```powershell
streamlit run Model_PichiaCLM/interfaces/streamlit_app.py --server.address 0.0.0.0 --server.port 8501
```

FastAPI mode:

1. Start the API.
2. Start Streamlit.
3. Select `FastAPI` in the sidebar and point it at `http://127.0.0.1:8000`.

## Configuration

The default model weights path is:

```text
Model_PichiaCLM/Training/PichiaData/2Target_AllData/Arch1-0404.weights.pt
```

Override it for API mode:

```powershell
$env:PICHIA_CLM_WEIGHTS = "C:\path\to\Arch1-0404.weights.pt"
$env:PICHIA_CLM_DEVICE = "cpu"
uvicorn Model_PichiaCLM.interfaces.api:app --host 127.0.0.1 --port 8000
```

By default, ambiguous amino acids such as `X`, `Z`, `B`, `U`, and `O` are rejected because they do not have a biological codon mask. Use `--allow-unknown` in the CLI or `"allow_unknown": true` in the API only when you intentionally want the notebook's unknown-token behavior.

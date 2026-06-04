# PichiaCLM-Torch

[中文说明](README.zh-CN.md) | English

PichiaCLM-Torch is a PyTorch implementation and deployable toolkit for codon optimization in *Pichia pastoris*.

The core workflow is:

```text
protein amino acid sequence -> optimized CDS/DNA sequence for Pichia expression
```

This project does not design new proteins. It converts an input amino acid sequence into a synonymous coding DNA sequence that better matches *Pichia pastoris* expression preferences, then provides sequence quality checks for cloning, synthesis, and review.

## Features

- Single amino acid sequence prediction.
- FASTA batch prediction for multiple proteins or variants.
- Streamlit web UI, FastAPI service, and CLI.
- Translation consistency checks for optimized CDS.
- GC and local GC analysis.
- CAI and codon usage comparison against both project training data and a public Kazusa *Pichia pastoris* codon table.
- Rare codon run, homopolymer, repeat, unwanted motif, and restriction enzyme site checks.
- External CDS quality control for sequences optimized outside this software.
- Signal peptide + mature protein comparison:
  - whole-sequence optimization
  - segmented optimization
- Conservative synonymous post-processing for selected sequence risks.

## Project Layout

```text
Model_PichiaCLM/
  core/
    predictor.py      # model loading and AA-to-CDS prediction
    analysis.py       # sequence quality analysis
    biology.py        # codon table, DNA normalization, translation checks
    fasta.py          # FASTA parsing and formatting
    restriction.py    # restriction enzyme site scanning
    fusion.py         # signal peptide + mature protein comparison
    postprocess.py    # conservative synonymous post-processing
  interfaces/
    cli.py            # command-line interface
    api.py            # FastAPI interface
    streamlit_app.py  # Streamlit web interface
```

The design keeps model inference in `core` and keeps CLI/API/Streamlit as thin interface layers.

## Installation

Core inference only:

```powershell
pip install -r requirements-core.txt
```

FastAPI:

```powershell
pip install -r requirements-api.txt
```

Streamlit:

```powershell
pip install -r requirements-streamlit.txt
```

Full deployment environment:

```powershell
pip install -r requirements-deploy.txt
```

## Quick Start

### CLI

Single prediction:

```powershell
python -m Model_PichiaCLM.interfaces.cli --aa MSTNPKPQR --json
```

Expected CDS for the bundled model:

```text
ATGTCCACAAATCCCAAACCACAGAGA
```

Batch FASTA prediction:

```powershell
python -m Model_PichiaCLM.interfaces.cli `
  --aa-fasta input.fasta `
  --analysis `
  --out-fasta output_cds.fasta `
  --out-csv report.csv
```

Analyze an externally optimized CDS without running model prediction:

```powershell
python -m Model_PichiaCLM.interfaces.cli `
  --cds ATGTCCACAAATCCCAAACCACAGAGA `
  --expected-aa MSTNPKPQR `
  --analysis
```

### FastAPI

```powershell
uvicorn Model_PichiaCLM.interfaces.api:app --host 0.0.0.0 --port 8000
```

Single prediction:

```powershell
Invoke-RestMethod -Method Post `
  -Uri http://127.0.0.1:8000/predict `
  -ContentType application/json `
  -Body '{"amino_acids":"MSTNPKPQR"}'
```

External CDS quality control:

```powershell
Invoke-RestMethod -Method Post `
  -Uri http://127.0.0.1:8000/analyze_cds `
  -ContentType application/json `
  -Body '{"cds":"ATGTCCACAAATCCCAAACCACAGAGA","expected_amino_acids":"MSTNPKPQR"}'
```

### Streamlit

```powershell
python -m streamlit run Model_PichiaCLM/interfaces/streamlit_app.py `
  --server.address 0.0.0.0 `
  --server.port 8501
```

Open:

```text
http://127.0.0.1:8501
```

The Streamlit UI includes:

- single prediction
- FASTA batch prediction
- external CDS quality control
- signal peptide fusion comparison
- conservative post-processing options

## Model Weights

Default weights path:

```text
Model_PichiaCLM/Training/PichiaData/2Target_AllData/Arch1-0404.weights.pt
```

The deployment code expects this file by default. You can override the path in CLI/API/Streamlit settings.

## Quality Checks

The quality report includes:

- CDS length and reading frame
- translation consistency against expected amino acids
- internal stop codons
- invalid DNA bases
- global GC and 30 bp sliding-window local GC
- CAI against training-data and public references
- codon usage statistics
- rare codon runs
- homopolymers and repeated sequence patterns
- default and custom restriction enzyme sites
- user-supplied unwanted motifs

Default GC thresholds:

```text
global GC: 35%-65%
local GC: 30 bp window, 25%-75%
```

Default restriction enzymes:

```text
EcoRI, XhoI, NotI, BamHI, HindIII, NdeI, NcoI, KpnI, XbaI, SpeI
```

## Notes

PichiaCLM output should be treated as a candidate optimized CDS. Before synthesis or experimental use, review translation consistency, cloning constraints, synthesis vendor rules, and project-specific biological requirements.

## Acknowledgments

This repository is based on the PyTorch port and deployment work for PichiaCLM-Torch. The original PichiaCLM concept and data-processing lineage come from prior PichiaCLM research code and datasets.

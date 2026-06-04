from __future__ import annotations

from dataclasses import asdict

import requests
import streamlit as st

from Model_PichiaCLM.core.config import DEFAULT_WEIGHTS_PATH
from Model_PichiaCLM.core.predictor import PichiaCLMPredictor


DEFAULT_SEQUENCE = "MSTNPKPQR"


@st.cache_resource(show_spinner="Loading PichiaCLM model...")
def load_predictor(weights_path: str, device: str | None) -> PichiaCLMPredictor:
    return PichiaCLMPredictor(weights_path=weights_path, device=device or None)


def predict_direct(amino_acids: str, allow_unknown: bool, weights_path: str, device: str | None):
    predictor = load_predictor(weights_path, device)
    return asdict(predictor.predict(amino_acids, allow_unknown=allow_unknown))


def predict_via_api(api_url: str, amino_acids: str, allow_unknown: bool):
    response = requests.post(
        f"{api_url.rstrip('/')}/predict",
        json={"amino_acids": amino_acids, "allow_unknown": allow_unknown},
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def main() -> None:
    st.set_page_config(page_title="PichiaCLM", page_icon="DNA", layout="wide")
    st.title("PichiaCLM AA-to-CDS")

    with st.sidebar:
        mode = st.radio("Prediction mode", ["Direct model", "FastAPI"], index=0)
        allow_unknown = st.checkbox("Allow ambiguous amino acids", value=False)
        if mode == "Direct model":
            weights_path = st.text_input("Weights path", value=str(DEFAULT_WEIGHTS_PATH))
            device = st.selectbox("Device", ["auto", "cpu", "cuda"], index=0)
            api_url = ""
        else:
            api_url = st.text_input("API URL", value="http://127.0.0.1:8000")
            weights_path = str(DEFAULT_WEIGHTS_PATH)
            device = "auto"

    amino_acids = st.text_area("Amino acid sequence", value=DEFAULT_SEQUENCE, height=160)

    if st.button("Predict", type="primary"):
        try:
            if mode == "Direct model":
                selected_device = None if device == "auto" else device
                result = predict_direct(amino_acids, allow_unknown, weights_path, selected_device)
            else:
                result = predict_via_api(api_url, amino_acids, allow_unknown)
        except Exception as exc:
            st.error(str(exc))
            return

        col_a, col_b, col_c = st.columns(3)
        col_a.metric("AA length", len(result["amino_acids"]))
        col_b.metric("CDS length", len(result["cds"]))
        col_c.metric("Device", result["device"])

        st.subheader("CDS")
        st.code(result["cds"], language="text")
        st.download_button(
            "Download CDS FASTA",
            data=f">PichiaCLM_prediction\n{result['cds']}\n",
            file_name="pichiaclm_prediction.fasta",
            mime="text/plain",
        )

        st.subheader("Codon IDs")
        st.json(result["codon_ids"])


if __name__ == "__main__":
    main()

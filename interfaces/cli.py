from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from Model_PichiaCLM.core.config import DEFAULT_WEIGHTS_PATH
from Model_PichiaCLM.core.predictor import PichiaCLMPredictor


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run PichiaCLM AA-to-CDS inference.")
    parser.add_argument("--aa", required=True, help="Amino acid sequence, e.g. MSTNPKPQR")
    parser.add_argument(
        "--weights",
        default=str(DEFAULT_WEIGHTS_PATH),
        help="Path to Arch1-0404.weights.pt",
    )
    parser.add_argument(
        "--device",
        default=None,
        help="Torch device override, e.g. cpu or cuda. Defaults to auto-detect.",
    )
    parser.add_argument(
        "--allow-unknown",
        action="store_true",
        help="Allow X/Z/B/U/O ambiguous amino acids using the notebook UNK token.",
    )
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    predictor = PichiaCLMPredictor(weights_path=args.weights, device=args.device)
    result = predictor.predict(args.aa, allow_unknown=args.allow_unknown)

    if args.json:
        print(json.dumps(asdict(result), ensure_ascii=False))
        return

    print(f"Input AA: {result.amino_acids}")
    print(f"Device: {result.device}")
    print(f"Codon IDs: {result.codon_ids}")
    print(f"CDS: {result.cds}")


if __name__ == "__main__":
    main()

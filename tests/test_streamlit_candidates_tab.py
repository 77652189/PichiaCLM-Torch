"""Headless end-to-end checks for the multi-candidate Streamlit tab.

These drive the real app through Streamlit's AppTest harness rather than a
browser, so widget state is committed deterministically. They exist because a
manual browser pass caught a NameError that no unit test could see: the
candidate result renderer referenced variables that were out of scope after a
refactor, which only surfaced once a full result was actually rendered.
"""
from __future__ import annotations

import unittest
from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1] / "Model_PichiaCLM" / "interfaces" / "streamlit_app.py"
RUN_TIMEOUT = 300

# 40 residues: long enough to fill a %MinMax window (18 codons) so the curve is
# actually plotted rather than short-circuited.
LONG_AA = "MSTNPKPQRSTNPKPQRTTNPKPQRSTNPKPQRTTNPKPQ"
# The same protein, coded with a different synonymous choice at every position:
# a valid already-aligned harmonization source per ADR-0008.
ALIGNED_SOURCE_CDS = "".join(
    {"M": "ATG", "S": "AGC", "T": "ACA", "N": "AAC", "P": "CCG", "K": "AAG", "Q": "CAG", "R": "CGC"}[aa]
    for aa in LONG_AA
)


def _candidates_tab_run(**widget_values) -> AppTest:
    app = AppTest.from_file(str(APP), default_timeout=RUN_TIMEOUT)
    app.run()

    app.text_area(key="candidate_aa").set_value(widget_values["amino_acids"])
    if widget_values.get("strategy_sampling"):
        app.radio(key="candidate_strategy").set_value("温度采样")
    if widget_values.get("harmonize"):
        app.checkbox(key="candidate_use_harmonization").set_value(True)
        app.run()
        app.text_area(key="candidate_source_cds").set_value(widget_values["source_cds"])

    app.button(key="candidate_predict").click()
    app.run()
    return app


class CandidatesTabRendersWithoutError(unittest.TestCase):
    def test_default_path_renders_a_full_result(self) -> None:
        app = _candidates_tab_run(amino_acids=LONG_AA)
        self.assertFalse(app.exception, [str(e) for e in app.exception])
        body = " ".join(str(e.value) for e in app.markdown) + " ".join(str(e.value) for e in app.caption)
        self.assertIn("排序依据", body)

    def test_harmonization_path_renders_a_full_result(self) -> None:
        app = _candidates_tab_run(
            amino_acids=LONG_AA,
            harmonize=True,
            source_cds=ALIGNED_SOURCE_CDS,
        )
        self.assertFalse(app.exception, [str(e) for e in app.exception])
        captions = " ".join(str(e.value) for e in app.caption)
        self.assertIn("harmonization", captions)

    def test_misaligned_source_cds_surfaces_an_actionable_error(self) -> None:
        """ADR-0008: a source CDS for a different protein must be refused in the
        UI with an explanation, not silently ignored or crash the page."""
        app = _candidates_tab_run(
            amino_acids=LONG_AA,
            harmonize=True,
            source_cds="ATGCCTCCTCCT",
        )
        self.assertFalse(app.exception, [str(e) for e in app.exception])
        errors = " ".join(str(e.value) for e in app.error)
        self.assertIn("does not encode the same protein", errors)


if __name__ == "__main__":
    unittest.main()

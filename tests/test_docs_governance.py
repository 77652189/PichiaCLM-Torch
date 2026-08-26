from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS = REPO_ROOT / "docs"
ADR_DIR = DOCS / "adr"
ACTIVE_ROOT_DOCS = {
    "ARCHITECTURE.md",
    "EXECUTION_PLAN.md",
    "HANDOFF.md",
    "REQUIREMENTS.md",
}
REQUIRED_HANDOFF_HEADINGS = {
    "## 当前切片",
    "## 下一步",
    "## 必读材料",
    "## 验证方式",
    "## 硬约束",
}


def test_docs_layout_is_the_reviewed_five_document_set():
    assert DOCS.is_dir()
    assert ADR_DIR.is_dir()
    assert {path.name for path in DOCS.glob("*.md")} == ACTIVE_ROOT_DOCS
    assert {path.name for path in ADR_DIR.glob("*.md")} == {
        "ADR-0001-qualified-candidate-acceptance.md",
        "ADR-0002-similarity-hard-threshold.md",
        "ADR-0003-min-max-host-only-profile.md",
        "ADR-0004-literature-informed-similarity-threshold.md",
        "ADR-0005-temperature-sampling-strategy-and-min-max-ranking.md",
        "ADR-0006-dynamic-source-reference-fetch.md",
        "ADR-0007-codon-axis-only-similarity-gate.md",
        "README.md",
    }


def test_adr_numbers_are_unique():
    """A repeated ADR number means two decisions claim the same identity, so
    every cross-reference to that number becomes ambiguous."""
    numbers = [path.name.split("-")[1] for path in ADR_DIR.glob("ADR-*.md")]
    assert numbers, "no ADR files found -- the glob or directory layout changed"
    assert len(numbers) == len(set(numbers)), f"duplicate ADR numbers: {sorted(numbers)}"


def test_handoff_preserves_required_sections_and_state_schema():
    handoff = (DOCS / "HANDOFF.md").read_text(encoding="utf-8")
    headings = {line.strip() for line in handoff.splitlines()}
    assert REQUIRED_HANDOFF_HEADINGS <= headings
    assert re.search(r"^current_slice: \S+$", handoff, re.MULTILINE)
    assert re.search(r"^slice_status: (planned|in_progress|done|blocked)$", handoff, re.MULTILINE)
    assert re.search(r"^authorization_status: (awaiting_user|authorized|paused)$", handoff, re.MULTILINE)
    assert re.search(r"^verification_status: (pending|passed|failed|not_run)$", handoff, re.MULTILINE)


def test_handoff_preserves_scientific_and_change_control_boundaries():
    handoff = (DOCS / "HANDOFF.md").read_text(encoding="utf-8")
    assert "不是表达产量、分泌效率或湿实验成功的预测" in handoff
    assert "不修改训练数据或模型权重" in handoff
    assert "CAI 不单独决定合格状态" in handoff
    assert "不改变远端可见性、不提交、不推送" in handoff


def test_adr_index_records_the_accepted_candidate_decision():
    index = (ADR_DIR / "README.md").read_text(encoding="utf-8")
    assert "[ADR-0001](ADR-0001-qualified-candidate-acceptance.md)" in index
    assert "accepted" in index

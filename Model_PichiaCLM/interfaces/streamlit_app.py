from __future__ import annotations

import csv
import importlib
import io
import json
from dataclasses import asdict
from typing import Any

import requests
import streamlit as st

from Model_PichiaCLM.core.analysis import analyze_cds, load_training_codon_reference
from Model_PichiaCLM.core.biology import normalize_dna
from Model_PichiaCLM.core.codon_editor import (
    build_codon_cells,
    codon_options_for_amino_acid,
    replace_selected_codons,
    search_codon_cells,
)
from Model_PichiaCLM.core.config import DEFAULT_WEIGHTS_PATH
from Model_PichiaCLM.core.fasta import FastaRecord, format_fasta, parse_fasta
from Model_PichiaCLM.core.fusion import compare_signal_fusion
from Model_PichiaCLM.core.postprocess import conservative_postprocess
from Model_PichiaCLM.core.predictor import PichiaCLMPredictor


DEFAULT_SEQUENCE = "MSTNPKPQR"
DIRECT_MODE_LABEL = "直接加载模型"
API_MODE_LABEL = "调用 FastAPI"


@st.cache_resource(show_spinner="正在加载 PichiaCLM 模型...")
def load_predictor(weights_path: str, device: str | None) -> PichiaCLMPredictor:
    predictor_cls = _current_predictor_class()
    return predictor_cls(weights_path=weights_path, device=device or None)


def _current_predictor_class() -> type[PichiaCLMPredictor]:
    module = importlib.import_module("Model_PichiaCLM.core.predictor")
    if not hasattr(module.PichiaCLMPredictor, "predict_candidates"):
        module = importlib.reload(module)
    return module.PichiaCLMPredictor


def load_predictor_with_candidates(weights_path: str, device: str | None) -> Any:
    predictor = load_predictor(weights_path, device)
    if hasattr(predictor, "predict_candidates"):
        return predictor
    load_predictor.clear()
    predictor = load_predictor(weights_path, device)
    if not hasattr(predictor, "predict_candidates"):
        raise RuntimeError("当前运行进程仍在使用旧版模型对象，请刷新页面或重启当前 Streamlit 预览进程。")
    return predictor


def parse_text_list(raw_text: str) -> list[str]:
    return [item.strip() for item in raw_text.replace(",", "\n").splitlines() if item.strip()]


def json_dumps(payload: object) -> str:
    return json.dumps(payload, ensure_ascii=False, indent=2)


def value_or_dash(value: object) -> object:
    return "-" if value is None else value


def gc_status_label(status: str) -> str:
    return {"ok": "正常", "low": "偏低", "high": "偏高"}.get(status, status)


def quality_status_label(status: str) -> str:
    return {"pass": "通过", "warning": "警告", "fail": "失败"}.get(status, status)


def candidate_source_label(source: str) -> str:
    return {
        "reference": "基准序列",
        "sample": "采样候选",
        "kazusa_constrained": "Kazusa 小幅替换",
        "kazusa_diverse": "Kazusa 多样性优化",
    }.get(source, source)


def quality_counts(analysis: dict[str, object]) -> tuple[int, int]:
    critical = 0
    if not analysis["valid_dna"]:
        critical += 1
    if not analysis["length_multiple_of_three"]:
        critical += 1
    if analysis["translation_matches_input"] is False:
        critical += 1
    critical += len(analysis["internal_stop_codons"])

    warnings = len(analysis["sequence_warnings"])
    warnings += 1 if analysis["gc_status"] != "ok" else 0
    warnings += len(analysis["local_gc_outliers"])
    warnings += len(analysis["restriction_sites"])
    warnings += len(analysis["motif_hits"])
    warnings += len(analysis["rare_codon_runs"])
    warnings += len(analysis["homopolymers"])
    warnings += len(analysis["tandem_repeats"])
    warnings += len(analysis["repeated_kmers"])
    return critical, warnings


def quality_status(analysis: dict[str, object]) -> str:
    critical, warnings = quality_counts(analysis)
    if critical:
        return "失败"
    if warnings:
        return "警告"
    return "通过"


def translation_status(analysis: dict[str, object]) -> str:
    if analysis["translation_matches_input"] is None:
        return "未提供参考 AA"
    return "通过" if analysis["translation_matches_input"] else "未通过"


def dataframe_or_success(title: str, rows: list[dict[str, object]], empty_text: str = "未发现") -> None:
    if rows:
        st.warning(f"{title}: 发现 {len(rows)} 项")
        st.dataframe(rows, use_container_width=True, hide_index=True)
    else:
        st.success(f"{title}: {empty_text}")


def records_to_csv(rows: list[dict[str, object]]) -> str:
    if not rows:
        return ""
    handle = io.StringIO()
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue()


def codon_display_value(cell: object, display_mode: str) -> str:
    return cell.rna_codon if display_mode == "RNA" else cell.dna_codon


def codon_option_display_value(row: dict[str, object], display_mode: str) -> str:
    return str(row["rna_codon"] if display_mode == "RNA" else row["codon"])


def render_codon_editor_styles() -> None:
    st.markdown(
        """
        <style>
        div[data-testid="stButton"] > button {
            min-height: 2.75rem;
            white-space: normal;
            border-radius: 0.55rem;
            font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
            line-height: 1.25;
            padding-top: 0.55rem;
            padding-bottom: 0.55rem;
        }
        div[data-testid="stButton"] > button[kind="primary"] {
            border: 2px solid #38bdf8;
            background: linear-gradient(180deg, #075985 0%, #0c4a6e 100%);
            box-shadow: 0 0 0 1px rgba(56, 189, 248, 0.35), 0 0 14px rgba(56, 189, 248, 0.25);
            color: #ffffff;
            font-weight: 700;
        }
        div[data-testid="stButton"] > button[kind="secondary"] {
            border: 1px solid rgba(148, 163, 184, 0.45);
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_codon_grid(
    cells: list[object],
    display_mode: str,
    matched_numbers: set[int],
    selected_numbers: set[int],
) -> None:
    st.caption("点击密码子块可选择或取消选择；蓝色高亮表示查找命中或已选择。")
    columns_per_row = 5
    version = st.session_state.get("codon_editor_selection_version", 0)
    for row_start in range(0, len(cells), columns_per_row):
        columns = st.columns(columns_per_row)
        for offset, cell in enumerate(cells[row_start : row_start + columns_per_row]):
            is_selected = cell.codon_number in selected_numbers
            is_matched = cell.codon_number in matched_numbers
            if is_selected:
                status = "已选"
            elif is_matched:
                status = "命中"
            elif not cell.replaceable:
                status = "不可替换"
            else:
                status = "可替换"
            label = (
                f"#{cell.codon_number}  {codon_display_value(cell, display_mode)}\n"
                f"AA {cell.amino_acid} | bp {cell.start}-{cell.end}\n"
                f"{status}"
            )
            clicked = columns[offset].button(
                label,
                key=f"codon_editor_cell_{version}_{cell.codon_number}",
                type="primary" if is_selected or is_matched else "secondary",
                disabled=not cell.replaceable,
                use_container_width=True,
            )
            if clicked:
                next_selected = set(selected_numbers)
                if is_selected:
                    next_selected.remove(cell.codon_number)
                else:
                    next_selected.add(cell.codon_number)
                st.session_state["codon_editor_selected_codons"] = sorted(next_selected)
                st.session_state["codon_editor_selection_version"] = version + 1
                st.rerun()


def codon_table_rows(
    cells: list[object],
    display_mode: str,
    matched_numbers: set[int],
    selected_numbers: set[int],
) -> list[dict[str, object]]:
    return [
        {
            "选择": cell.codon_number in selected_numbers,
            "命中": "是" if cell.codon_number in matched_numbers else "",
            "密码子编号": cell.codon_number,
            "bp 位置": f"{cell.start}-{cell.end}",
            "显示密码子": codon_display_value(cell, display_mode),
            "DNA": cell.dna_codon,
            "RNA": cell.rna_codon,
            "氨基酸": cell.amino_acid,
            "可替换": "是" if cell.replaceable else "",
            "训练频率": cell.training_fraction,
            "公开表频率": cell.public_fraction,
        }
        for cell in cells
    ]


def replacement_rows(replacements: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        {
            "密码子编号": row["codon_number"],
            "bp 位置": f"{row['start']}-{row['end']}",
            "氨基酸": row["amino_acid"],
            "原密码子": row["old_codon"],
            "新密码子": row["new_codon"],
        }
        for row in replacements
    ]


def predict_direct(
    amino_acids: str,
    allow_unknown: bool,
    weights_path: str,
    device: str | None,
    motifs: list[str],
    custom_sites: list[str],
    do_postprocess: bool,
) -> dict[str, object]:
    predictor = load_predictor(weights_path, device)
    result = predictor.predict(amino_acids, allow_unknown=allow_unknown)
    payload = asdict(result)
    payload["analysis"] = asdict(
        analyze_cds(
            result.cds,
            amino_acids=result.amino_acids,
            motifs=motifs,
            custom_restriction_sites=custom_sites,
        )
    )
    if do_postprocess:
        training_reference, _ = load_training_codon_reference()
        payload["postprocess"] = asdict(
            conservative_postprocess(
                result.cds,
                result.amino_acids,
                reference_fractions=training_reference,
                forbidden_motifs=motifs,
                custom_restriction_sites=custom_sites,
            )
        )
    return payload


def predict_via_api(
    api_url: str,
    amino_acids: str,
    allow_unknown: bool,
    motifs: list[str],
    custom_sites: list[str],
    do_postprocess: bool,
) -> dict[str, object]:
    response = requests.post(
        f"{api_url.rstrip('/')}/predict",
        json={
            "amino_acids": amino_acids,
            "allow_unknown": allow_unknown,
            "include_analysis": True,
            "unwanted_motifs": motifs,
            "custom_restriction_sites": custom_sites,
            "postprocess": do_postprocess,
        },
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def predict_candidates_direct(
    amino_acids: str,
    allow_unknown: bool,
    weights_path: str,
    device: str | None,
    motifs: list[str],
    custom_sites: list[str],
    num_candidates: int,
    subset_size: int,
    temperature: float,
    seed: int | None,
) -> dict[str, object]:
    predictor = load_predictor_with_candidates(weights_path, device)
    return asdict(
        predictor.predict_candidates(
            amino_acids,
            allow_unknown=allow_unknown,
            num_candidates=num_candidates,
            subset_size=subset_size,
            temperature=temperature,
            seed=seed,
            motifs=motifs,
            custom_restriction_sites=custom_sites,
        )
    )


def predict_candidates_via_api(
    api_url: str,
    amino_acids: str,
    allow_unknown: bool,
    motifs: list[str],
    custom_sites: list[str],
    num_candidates: int,
    subset_size: int,
    temperature: float,
    seed: int | None,
) -> dict[str, object]:
    response = requests.post(
        f"{api_url.rstrip('/')}/predict_candidates",
        json={
            "amino_acids": amino_acids,
            "num_candidates": num_candidates,
            "subset_size": subset_size,
            "temperature": temperature,
            "seed": seed,
            "allow_unknown": allow_unknown,
            "unwanted_motifs": motifs,
            "custom_restriction_sites": custom_sites,
        },
        timeout=180,
    )
    response.raise_for_status()
    return response.json()


def analyze_external_cds_direct(
    cds: str,
    expected_amino_acids: str | None,
    motifs: list[str],
    custom_sites: list[str],
) -> dict[str, object]:
    analysis = analyze_cds(
        cds,
        amino_acids=expected_amino_acids,
        motifs=motifs,
        custom_restriction_sites=custom_sites,
    )
    return {
        "cds": normalize_dna(cds),
        "expected_amino_acids": expected_amino_acids,
        "translated_amino_acids": analysis.translated_amino_acids,
        "analysis": asdict(analysis),
    }


def analyze_external_cds_via_api(
    api_url: str,
    cds: str,
    expected_amino_acids: str | None,
    motifs: list[str],
    custom_sites: list[str],
) -> dict[str, object]:
    response = requests.post(
        f"{api_url.rstrip('/')}/analyze_cds",
        json={
            "cds": cds,
            "expected_amino_acids": expected_amino_acids,
            "unwanted_motifs": motifs,
            "custom_restriction_sites": custom_sites,
        },
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def run_prediction(
    mode: str,
    amino_acids: str,
    allow_unknown: bool,
    weights_path: str,
    device: str | None,
    api_url: str,
    motifs: list[str],
    custom_sites: list[str],
    do_postprocess: bool,
) -> dict[str, object]:
    if mode == DIRECT_MODE_LABEL:
        return predict_direct(amino_acids, allow_unknown, weights_path, device, motifs, custom_sites, do_postprocess)
    return predict_via_api(api_url, amino_acids, allow_unknown, motifs, custom_sites, do_postprocess)


def run_candidate_prediction(
    mode: str,
    amino_acids: str,
    allow_unknown: bool,
    weights_path: str,
    device: str | None,
    api_url: str,
    motifs: list[str],
    custom_sites: list[str],
    num_candidates: int,
    subset_size: int,
    temperature: float,
    seed: int | None,
) -> dict[str, object]:
    if mode == DIRECT_MODE_LABEL:
        return predict_candidates_direct(
            amino_acids,
            allow_unknown,
            weights_path,
            device,
            motifs,
            custom_sites,
            num_candidates,
            subset_size,
            temperature,
            seed,
        )
    return predict_candidates_via_api(
        api_url,
        amino_acids,
        allow_unknown,
        motifs,
        custom_sites,
        num_candidates,
        subset_size,
        temperature,
        seed,
    )


def run_cds_analysis(
    mode: str,
    cds: str,
    expected_amino_acids: str | None,
    api_url: str,
    motifs: list[str],
    custom_sites: list[str],
) -> dict[str, object]:
    if mode == DIRECT_MODE_LABEL:
        return analyze_external_cds_direct(cds, expected_amino_acids, motifs, custom_sites)
    return analyze_external_cds_via_api(api_url, cds, expected_amino_acids, motifs, custom_sites)


def render_quality_overview(analysis: dict[str, object]) -> None:
    critical, warnings = quality_counts(analysis)
    status = quality_status(analysis)
    if status == "通过":
        st.success("总体结论: 通过当前质量检查")
    elif status == "警告":
        st.warning("总体结论: 存在需要复核的风险项")
    else:
        st.error("总体结论: 存在基础正确性问题")

    col_a, col_b, col_c = st.columns(3)
    col_a.metric("结论", status)
    col_b.metric("基础问题", critical)
    col_c.metric("风险警告", warnings)


def render_quality_report(analysis: dict[str, object]) -> None:
    st.subheader("质量警告")
    with st.expander("基础正确性检查", expanded=True):
        for warning in analysis["sequence_warnings"]:
            st.warning(warning)
        st.success("DNA 字母检查: 通过" if analysis["valid_dna"] else "DNA 字母检查: 未通过")
        st.success("阅读框检查: 通过" if analysis["length_multiple_of_three"] else "阅读框检查: 未通过")
        if analysis["translation_matches_input"] is None:
            st.info("翻译一致性: 未提供参考 AA")
        else:
            st.success("翻译一致性: 通过" if analysis["translation_matches_input"] else "翻译一致性: 未通过")
        dataframe_or_success(
            "内部终止密码子",
            [{"密码子编号": codon_number} for codon_number in analysis["internal_stop_codons"]],
        )

    with st.expander("GC 含量"):
        st.write(f"全局 GC 含量: {analysis['gc_percent']}% ({gc_status_label(analysis['gc_status'])})")
        dataframe_or_success(
            "局部 GC 异常区域",
            [
                {"起始 bp": row["start"], "结束 bp": row["end"], "GC (%)": row["gc_percent"]}
                for row in analysis["local_gc_outliers"][:50]
            ],
        )

    with st.expander("酶切位点与 motif"):
        dataframe_or_success(
            "限制性酶切位点",
            [
                {
                    "酶名": row["name"],
                    "识别序列": row["sequence"],
                    "起始 bp": row["start"],
                    "结束 bp": row["end"],
                }
                for row in analysis["restriction_sites"]
            ],
        )
        dataframe_or_success(
            "不想要的 motif",
            [
                {"Motif 序列": row["motif"], "起始 bp": row["start"], "结束 bp": row["end"]}
                for row in analysis["motif_hits"]
            ],
        )

    with st.expander("密码子与重复序列"):
        dataframe_or_success(
            "连续稀有密码子",
            [
                {
                    "参考表": "训练数据" if row["reference"] == "training" else "公开表",
                    "起始密码子": row["start_codon"],
                    "结束密码子": row["end_codon"],
                    "密码子": " ".join(row["codons"]),
                }
                for row in analysis["rare_codon_runs"][:50]
            ],
        )
        dataframe_or_success(
            "长同聚碱基",
            [
                {"碱基": row["base"], "起始 bp": row["start"], "结束 bp": row["end"], "长度": row["length"]}
                for row in analysis["homopolymers"][:50]
            ],
        )
        dataframe_or_success(
            "串联重复",
            [
                {
                    "重复单元": row["sequence"],
                    "起始 bp": row["start"],
                    "结束 bp": row["end"],
                    "重复次数": row["copies"],
                }
                for row in analysis["tandem_repeats"][:50]
            ],
        )
        dataframe_or_success(
            "重复 12 bp 片段",
            [
                {
                    "重复片段": row["sequence"],
                    "出现次数": row["count"],
                    "前几个位置": ", ".join(str(position) for position in row["positions"]),
                }
                for row in analysis["repeated_kmers"][:50]
            ],
        )


def render_codon_usage(analysis: dict[str, object]) -> None:
    with st.expander("密码子使用对比"):
        used_rows = [
            {
                "密码子": row["codon"],
                "氨基酸": row["amino_acid"],
                "出现次数": row["count"],
                "本序列比例": row["sequence_fraction"],
                "训练数据比例": row["training_fraction"],
                "公开表比例": row["public_fraction"],
            }
            for row in analysis["codon_usage"]
            if row["count"] > 0
        ]
        st.dataframe(used_rows, use_container_width=True, hide_index=True)


def render_postprocess(result: dict[str, object], key_prefix: str) -> None:
    postprocess = result.get("postprocess")
    if not postprocess:
        return
    st.subheader("后处理建议")
    st.metric("同义替换次数", len(postprocess["replacements"]))
    st.success("翻译保持一致" if postprocess["translation_preserved"] else "翻译未保持一致")
    if postprocess["optimized_cds"] != postprocess["original_cds"]:
        st.text_area("后处理后的 CDS", value=postprocess["optimized_cds"], height=120, key=f"{key_prefix}_post_cds")
    dataframe_or_success(
        "替换记录",
        [
            {
                "密码子编号": row["codon_number"],
                "位置": f"{row['start']}-{row['end']}",
                "氨基酸": row["amino_acid"],
                "原密码子": row["old_codon"],
                "新密码子": row["new_codon"],
                "原因": row["reason"],
            }
            for row in postprocess["replacements"]
        ],
    )
    dataframe_or_success("仍未解决的问题", [{"问题": item} for item in postprocess["remaining_issues"]])


def render_prediction_result(result: dict[str, object], title: str = "预测结果", key_prefix: str = "result") -> None:
    analysis = result["analysis"]
    st.subheader(title)
    col_a, col_b, col_c = st.columns(3)
    col_a.metric("氨基酸长度", len(result["amino_acids"]))
    col_b.metric("CDS 长度", len(result["cds"]))
    col_c.metric("运行设备", result["device"])

    st.text_area("优化后的 CDS", value=result["cds"], height=120, key=f"{key_prefix}_cds_text")
    st.download_button(
        "下载 CDS FASTA",
        data=f">PichiaCLM_prediction\n{result['cds']}\n",
        file_name="pichiaclm_prediction.fasta",
        mime="text/plain",
        key=f"{key_prefix}_cds_download",
    )

    metric_a, metric_b, metric_c, metric_d = st.columns(4)
    metric_a.metric("GC 含量", f"{analysis['gc_percent']}%", gc_status_label(analysis["gc_status"]))
    metric_b.metric("局部 GC 警告", len(analysis["local_gc_outliers"]))
    metric_c.metric("CAI 训练数据", value_or_dash(analysis["cai"]["training"]))
    metric_d.metric("CAI 公开表", value_or_dash(analysis["cai"]["public"]))

    check_a, check_b, check_c, check_d = st.columns(4)
    check_a.metric("翻译一致性", translation_status(analysis))
    check_b.metric("酶切位点", len(analysis["restriction_sites"]))
    check_c.metric("Motif 命中", len(analysis["motif_hits"]))
    check_d.metric("连续稀有密码子", len(analysis["rare_codon_runs"]))

    st.caption("GC 阈值: 全局 35%-65%；局部 30 bp 窗口 25%-75%。公开 CAI 参考表: Kazusa Pichia pastoris taxon 4922。")
    render_quality_overview(analysis)
    render_quality_report(analysis)
    render_codon_usage(analysis)
    render_postprocess(result, key_prefix=key_prefix)

    st.download_button(
        "下载分析报告 JSON",
        data=json_dumps(result),
        file_name="pichiaclm_analysis.json",
        mime="application/json",
        key=f"{key_prefix}_json_download",
    )


def render_cds_analysis_result(result: dict[str, object], title: str = "CDS 质检结果", key_prefix: str = "cds_qc") -> None:
    analysis = result["analysis"]
    st.subheader(title)
    col_a, col_b, col_c = st.columns(3)
    col_a.metric("CDS 长度", analysis["cds_length"])
    col_b.metric("密码子数量", analysis["codon_count"])
    col_c.metric("翻译一致性", translation_status(analysis))

    st.text_area("待质检 CDS", value=result["cds"], height=120, key=f"{key_prefix}_cds_text")
    st.text_area("翻译得到的 AA", value=result["translated_amino_acids"], height=100, key=f"{key_prefix}_translated_aa")

    metric_a, metric_b, metric_c, metric_d = st.columns(4)
    metric_a.metric("GC 含量", f"{analysis['gc_percent']}%", gc_status_label(analysis["gc_status"]))
    metric_b.metric("局部 GC 警告", len(analysis["local_gc_outliers"]))
    metric_c.metric("CAI 训练数据", value_or_dash(analysis["cai"]["training"]))
    metric_d.metric("CAI 公开表", value_or_dash(analysis["cai"]["public"]))

    check_a, check_b, check_c = st.columns(3)
    check_a.metric("酶切位点", len(analysis["restriction_sites"]))
    check_b.metric("Motif 命中", len(analysis["motif_hits"]))
    check_c.metric("非法碱基", len(analysis["invalid_bases"]))

    render_quality_overview(analysis)
    render_quality_report(analysis)
    render_codon_usage(analysis)
    st.download_button(
        "下载 CDS 质检报告 JSON",
        data=json_dumps(result),
        file_name="pichiaclm_cds_qc.json",
        mime="application/json",
        key=f"{key_prefix}_json_download",
    )


def render_candidate_set_result(result: dict[str, object]) -> None:
    diversity = result["pairwise_diversity"]
    recommended_subset = result.get("recommended_subset")
    recommended_ranks = set(recommended_subset["selected_ranks"]) if recommended_subset else set()
    st.subheader("候选对比")
    col_a, col_b, col_c, col_d = st.columns(4)
    col_a.metric("生成候选", f"{result['generated_candidates']} / {result['requested_candidates']}")
    col_b.metric("采样次数", result["attempts"])
    col_c.metric("平均密码子差异", value_or_dash(diversity["mean_codon_difference_percent"]))
    col_d.metric("最小密码子差异", value_or_dash(diversity["min_codon_difference_percent"]))
    if result.get("note"):
        st.info(result["note"])

    rows = []
    for candidate in result["candidates"]:
        analysis = candidate["analysis"]
        quality = candidate["quality"]
        difference = candidate["difference_from_reference"]
        preference = candidate["codon_preference"]
        rows.append(
            {
                "排名": candidate["rank"],
                "推荐子集": "是" if candidate["rank"] in recommended_ranks else "",
                "来源": candidate_source_label(candidate["source"]),
                "结论": quality_status_label(quality["status"]),
                "基础问题": quality["critical_issues"],
                "风险警告": quality["warnings"],
                "GC%": analysis["gc_percent"],
                "CAI 训练数据": analysis["cai"]["training"],
                "CAI 公开表": analysis["cai"]["public"],
                "bp 差异%": difference["bp_difference_percent"],
                "密码子差异%": difference["codon_difference_percent"],
                "Kazusa 最优%": preference["top_preferred_percent"],
                "Kazusa 次优%": preference["second_preferred_percent"],
                "Kazusa 最低频%": preference["lowest_preferred_percent"],
                "酶切位点": len(analysis["restriction_sites"]),
                "Motif 命中": len(analysis["motif_hits"]),
                "局部 GC 警告": len(analysis["local_gc_outliers"]),
            }
        )
    st.dataframe(rows, use_container_width=True, hide_index=True)

    pairwise_rows = result.get("pairwise_similarities", [])
    if pairwise_rows:
        similarity_values = [row["codon_similarity_percent"] for row in pairwise_rows]
        st.subheader("候选间相似度")
        sim_a, sim_b, sim_c = st.columns(3)
        sim_a.metric("平均密码子相似度", value_or_dash(round(sum(similarity_values) / len(similarity_values), 2)))
        sim_b.metric("最高密码子相似度", value_or_dash(max(similarity_values)))
        sim_c.metric("最低密码子相似度", value_or_dash(min(similarity_values)))
        st.caption("相似度越低，说明两条 CDS 在同义密码子选择上越不一样；所有序列仍应翻译回同一个蛋白。")
        similarity_table = [
            {
                "候选 A": row["left_rank"],
                "候选 B": row["right_rank"],
                "密码子相似度%": row["codon_similarity_percent"],
                "密码子差异%": row["codon_difference_percent"],
                "bp 相似度%": row["bp_similarity_percent"],
                "bp 差异%": row["bp_difference_percent"],
            }
            for row in pairwise_rows
        ]
        st.dataframe(similarity_table, use_container_width=True, hide_index=True)

    if recommended_subset:
        st.subheader("推荐低相似度子集")
        subset_a, subset_b, subset_c = st.columns(3)
        subset_a.metric("推荐条数", f"{recommended_subset['selected_size']} / {recommended_subset['requested_size']}")
        subset_b.metric("平均密码子相似度", value_or_dash(recommended_subset["mean_codon_similarity_percent"]))
        subset_c.metric("最高密码子相似度", value_or_dash(recommended_subset["max_codon_similarity_percent"]))
        st.write("推荐候选排名：" + "、".join(str(rank) for rank in recommended_subset["selected_ranks"]))
        st.caption("这组序列是在已生成候选中挑出的低相似度组合；所有序列仍保留各自的翻译一致性和质量检查结果。")
        subset_pair_rows = [
            {
                "候选 A": row["left_rank"],
                "候选 B": row["right_rank"],
                "密码子相似度%": row["codon_similarity_percent"],
                "密码子差异%": row["codon_difference_percent"],
                "bp 相似度%": row["bp_similarity_percent"],
                "bp 差异%": row["bp_difference_percent"],
            }
            for row in pairwise_rows
            if row["left_rank"] in recommended_ranks and row["right_rank"] in recommended_ranks
        ]
        if subset_pair_rows:
            st.dataframe(subset_pair_rows, use_container_width=True, hide_index=True)
        selected_candidates = [
            candidate for candidate in result["candidates"]
            if candidate["rank"] in recommended_ranks
        ]
        selected_fasta_records = [
            FastaRecord(
                id=f"candidate_{candidate['rank']}_{candidate['source']}",
                description="PichiaCLM low-similarity CDS subset",
                sequence=candidate["cds"],
            )
            for candidate in selected_candidates
        ]
        selected_rows = [row for row in rows if row["排名"] in recommended_ranks]
        subset_fasta_col, subset_csv_col = st.columns(2)
        subset_fasta_col.download_button(
            "下载推荐子集 FASTA",
            data=format_fasta(selected_fasta_records),
            file_name="pichiaclm_low_similarity_subset.fasta",
            mime="text/plain",
        )
        subset_csv_col.download_button(
            "下载推荐子集 CSV",
            data=records_to_csv(selected_rows),
            file_name="pichiaclm_low_similarity_subset.csv",
            mime="text/csv",
        )

    fasta_records = [
        FastaRecord(
            id=f"candidate_{candidate['rank']}_{candidate['source']}",
            description="PichiaCLM CDS candidate",
            sequence=candidate["cds"],
        )
        for candidate in result["candidates"]
    ]
    col_fasta, col_csv = st.columns(2)
    col_fasta.download_button(
        "下载候选 CDS FASTA",
        data=format_fasta(fasta_records),
        file_name="pichiaclm_candidates.fasta",
        mime="text/plain",
    )
    col_csv.download_button(
        "下载候选对比 CSV",
        data=records_to_csv(rows),
        file_name="pichiaclm_candidates.csv",
        mime="text/csv",
    )

    with st.expander("逐条候选详情", expanded=False):
        for candidate in result["candidates"]:
            analysis = candidate["analysis"]
            difference = candidate["difference_from_reference"]
            preference = candidate["codon_preference"]
            st.markdown(f"### 候选 {candidate['rank']} - {candidate_source_label(candidate['source'])}")
            st.text_area(
                "CDS",
                value=candidate["cds"],
                height=110,
                key=f"candidate_{candidate['rank']}_cds",
            )
            st.write(
                f"GC {analysis['gc_percent']}%；CAI 训练数据 {value_or_dash(analysis['cai']['training'])}；"
                f"与基准差异 {difference['bp_differences']} bp / {difference['codon_differences']} 个密码子。"
            )
            st.write(
                f"Kazusa 偏好排名：最优 {preference['top_preferred_percent']}%，"
                f"次优 {preference['second_preferred_percent']}%，"
                f"最低频 {preference['lowest_preferred_percent']}%。"
            )
            render_quality_overview(analysis)

    st.download_button(
        "下载完整候选报告 JSON",
        data=json_dumps(result),
        file_name="pichiaclm_candidates.json",
        mime="application/json",
    )


def sidebar_settings() -> dict[str, object]:
    with st.sidebar:
        mode = st.radio("预测模式", [DIRECT_MODE_LABEL, API_MODE_LABEL], index=0)
        if mode == DIRECT_MODE_LABEL:
            weights_path = st.text_input("模型权重路径", value=str(DEFAULT_WEIGHTS_PATH))
            device = st.selectbox("运行设备", ["auto", "cpu", "cuda"], index=0)
            api_url = ""
        else:
            api_url = st.text_input("API 服务地址", value="http://127.0.0.1:8000")
            weights_path = str(DEFAULT_WEIGHTS_PATH)
            device = "auto"
        with st.expander("高级预测设置"):
            allow_unknown = st.checkbox("允许模糊氨基酸", value=False)
            do_postprocess = st.checkbox("启用保守后处理", value=False)
        with st.expander("高级质检设置"):
            motifs = parse_text_list(st.text_area("不想要的 motif", value="", height=90))
            custom_sites = parse_text_list(st.text_area("自定义酶切位点", value="", height=90))
    return {
        "mode": mode,
        "allow_unknown": allow_unknown,
        "do_postprocess": do_postprocess,
        "motifs": motifs,
        "custom_sites": custom_sites,
        "weights_path": weights_path,
        "device": None if device == "auto" else device,
        "api_url": api_url,
    }


def render_single_tab(settings: dict[str, object]) -> None:
    amino_acids = st.text_area("氨基酸序列", value=DEFAULT_SEQUENCE, height=160)
    if st.button("开始预测", type="primary", key="single_predict"):
        try:
            result = run_prediction(amino_acids=amino_acids, **settings)
        except Exception as exc:
            st.error(str(exc))
            return
        render_prediction_result(result)


def render_candidates_tab(settings: dict[str, object]) -> None:
    amino_acids = st.text_area("氨基酸序列", value=DEFAULT_SEQUENCE, height=140, key="candidate_aa")
    col_a, col_b, col_c, col_d = st.columns(4)
    num_candidates = int(col_a.number_input("候选数量", min_value=2, max_value=50, value=10, step=1))
    subset_size = int(
        col_b.number_input(
            "推荐子集大小",
            min_value=1,
            max_value=num_candidates,
            value=min(5, num_candidates),
            step=1,
        )
    )
    temperature = float(col_c.slider("采样温度", min_value=0.1, max_value=2.0, value=0.8, step=0.1))
    use_seed = col_d.checkbox("固定随机种子", value=True)
    seed = int(st.number_input("随机种子", min_value=0, value=42, step=1)) if use_seed else None
    st.caption(
        "候选序列都必须翻译回同一个 AA。当前策略从基准 CDS 出发做 10%-20% 以内的小幅同义替换，"
        "优先使用 Kazusa 次优/中频密码子，并尽量不增加最低频密码子。"
    )
    if settings["do_postprocess"]:
        st.info("多候选 CDS 默认不自动执行保守后处理；建议先选择候选，再进入单条或二次 CDS 质检流程复核。")

    if st.button("生成候选 CDS", type="primary", key="candidate_predict"):
        try:
            result = run_candidate_prediction(
                mode=settings["mode"],
                amino_acids=amino_acids,
                allow_unknown=settings["allow_unknown"],
                weights_path=settings["weights_path"],
                device=settings["device"],
                api_url=settings["api_url"],
                motifs=settings["motifs"],
                custom_sites=settings["custom_sites"],
                num_candidates=num_candidates,
                subset_size=subset_size,
                temperature=temperature,
                seed=seed,
            )
        except Exception as exc:
            st.error(str(exc))
            return
        render_candidate_set_result(result)


def render_codon_editor_tab(settings: dict[str, object]) -> None:
    render_codon_editor_styles()
    default_cds = "ATGTCCACAAATCCCAAACCACAGAGA"
    st.session_state.setdefault("codon_editor_active_cds", default_cds)
    st.session_state.setdefault("codon_editor_selected_codons", [])
    st.session_state.setdefault("codon_editor_selection_version", 0)
    st.session_state.setdefault("codon_editor_search_query", "")

    input_col, option_col = st.columns([2, 1])
    with input_col:
        raw_cds = st.text_area("粘贴 CDS", value=default_cds, height=140, key="codon_editor_raw_cds")
        expected_aa = st.text_area("可选：期望翻译得到的 AA", value="", height=95, key="codon_editor_expected_aa")
    with option_col:
        show_rna = st.toggle("显示 RNA 写法 U", value=False)
        display_mode = "RNA" if show_rna else "DNA"
        st.caption("编辑器内部始终按 DNA/CDS 处理；RNA 显示只用于对照标准密码子表。")
        if st.button("载入到编辑器", type="primary", key="codon_editor_load"):
            try:
                build_codon_cells(raw_cds)
            except ValueError as exc:
                st.error(str(exc))
            else:
                st.session_state["codon_editor_active_cds"] = normalize_dna(raw_cds)
                st.session_state["codon_editor_selected_codons"] = []
                st.session_state["codon_editor_search_query"] = ""
                st.session_state["codon_editor_last_result"] = None
                st.session_state["codon_editor_selection_version"] += 1
                st.rerun()

    active_cds = st.session_state["codon_editor_active_cds"]
    try:
        cells = build_codon_cells(active_cds)
    except ValueError as exc:
        st.error(str(exc))
        return

    st.subheader("当前编辑中的 CDS")
    st.text_area(
        "当前 CDS",
        value=active_cds,
        height=110,
        disabled=True,
        key=f"codon_editor_active_cds_display_{st.session_state['codon_editor_selection_version']}",
    )

    query_col, apply_col = st.columns([3, 1])
    query = query_col.text_input("查找密码子或氨基酸", value="", placeholder="例如 CCA / UUC / P / Phe")
    apply_col.write("")
    apply_col.write("")
    if apply_col.button("应用查找", type="primary", key="codon_editor_apply_search"):
        st.session_state["codon_editor_search_query"] = query.strip()
        st.session_state["codon_editor_selection_version"] += 1
        st.rerun()
    typed_query = query.strip()
    active_query = typed_query or st.session_state.get("codon_editor_search_query", "")
    matched_numbers = set(search_codon_cells(cells, active_query))
    selected_numbers = set(st.session_state.get("codon_editor_selected_codons", []))

    metric_a, metric_b, metric_c, metric_d = st.columns(4)
    metric_a.metric("密码子数量", len(cells))
    metric_b.metric("命中数量", len(matched_numbers))
    metric_c.metric("已选择", len(selected_numbers))
    metric_d.metric("显示模式", display_mode)

    if active_query:
        if matched_numbers:
            matched_preview = ", ".join(
                f"#{cell.codon_number} {codon_display_value(cell, display_mode)}"
                for cell in cells
                if cell.codon_number in matched_numbers
            )
            st.success(f"当前查找：{active_query}，命中 {len(matched_numbers)} 个：{matched_preview}")
        else:
            st.warning(f"当前查找：{active_query}，没有命中。")
    elif query.strip():
        st.info("输入查找内容后，点击“应用查找”刷新高亮。")

    select_col, clear_col = st.columns(2)
    if select_col.button("全选当前命中", disabled=not matched_numbers, key="codon_editor_select_matches"):
        st.session_state["codon_editor_selected_codons"] = sorted(selected_numbers | matched_numbers)
        st.session_state["codon_editor_selection_version"] += 1
        st.rerun()
    if clear_col.button("清空选择", disabled=not selected_numbers, key="codon_editor_clear_selection"):
        st.session_state["codon_editor_selected_codons"] = []
        st.session_state["codon_editor_selection_version"] += 1
        st.rerun()

    render_codon_grid(cells, display_mode, matched_numbers, selected_numbers)

    rows = codon_table_rows(cells, display_mode, matched_numbers, selected_numbers)
    with st.expander("详细表格：逐个勾选密码子", expanded=False):
        edited_rows = st.data_editor(
            rows,
            use_container_width=True,
            hide_index=True,
            key=f"codon_editor_table_{st.session_state['codon_editor_selection_version']}_{len(cells)}",
            disabled=[column for column in rows[0] if column != "选择"] if rows else True,
            column_config={
                "选择": st.column_config.CheckboxColumn("选择"),
            },
        )
    selected_numbers = {
        int(row["密码子编号"])
        for row in edited_rows
        if row.get("选择")
    }
    st.session_state["codon_editor_selected_codons"] = sorted(selected_numbers)

    selected_cells = [
        cell for cell in cells
        if cell.codon_number in selected_numbers
    ]
    selected_amino_acids = sorted({cell.amino_acid for cell in selected_cells})

    st.subheader("同义密码子替换")
    if not selected_cells:
        st.info("请先查找并选择至少一个密码子。")
    elif len(selected_amino_acids) > 1:
        st.warning("当前选择包含多个氨基酸，请一次只替换同一种氨基酸对应的密码子。")
    else:
        amino_acid = selected_amino_acids[0]
        options = codon_options_for_amino_acid(amino_acid)
        if amino_acid == "*" or len(options) <= 1:
            st.warning(f"{amino_acid} 没有可替换的同义密码子。")
        else:
            option_labels = {
                (
                    f"{codon_option_display_value(row, display_mode)} "
                    f"| 训练频率 {value_or_dash(row['training_fraction'])} "
                    f"| 公开表频率 {value_or_dash(row['public_fraction'])}"
                ): row["codon"]
                for row in options
            }
            target_label = st.selectbox(
                "目标同义密码子",
                list(option_labels),
                key=f"codon_editor_target_codon_{display_mode}_{amino_acid}",
            )
            target_display = target_label.split(" | ", maxsplit=1)[0]
            st.caption(f"将 {len(selected_cells)} 个 {amino_acid} 密码子替换为 {target_display}。")
            if st.button("确认替换选中密码子", type="primary", key="codon_editor_replace"):
                try:
                    edit_result = replace_selected_codons(
                        active_cds,
                        selected_numbers,
                        option_labels[target_label],
                        expected_amino_acids=expected_aa or None,
                    )
                    analysis = analyze_cds(
                        edit_result.edited_cds,
                        amino_acids=edit_result.expected_amino_acids,
                        motifs=settings["motifs"],
                        custom_restriction_sites=settings["custom_sites"],
                    )
                except ValueError as exc:
                    st.error(str(exc))
                else:
                    st.session_state["codon_editor_active_cds"] = edit_result.edited_cds
                    st.session_state["codon_editor_selected_codons"] = []
                    st.session_state["codon_editor_selection_version"] += 1
                    st.session_state["codon_editor_last_result"] = {
                        "edit": asdict(edit_result),
                        "analysis": asdict(analysis),
                    }
                    st.rerun()

    last_result = st.session_state.get("codon_editor_last_result")
    if last_result:
        edit = last_result["edit"]
        st.subheader("替换后结果")
        st.success("翻译保持一致" if edit["translation_preserved"] else "翻译未保持一致")
        rows = replacement_rows(edit["replacements"])
        dataframe_or_success("替换记录", rows, empty_text="没有实际替换")
        download_col_a, download_col_b = st.columns(2)
        download_col_a.download_button(
            "下载替换后 CDS FASTA",
            data=f">PichiaCLM_codon_edited\n{edit['edited_cds']}\n",
            file_name="pichiaclm_codon_edited.fasta",
            mime="text/plain",
            key="codon_editor_fasta_download",
        )
        download_col_b.download_button(
            "下载替换记录 CSV",
            data=records_to_csv(rows),
            file_name="pichiaclm_codon_replacements.csv",
            mime="text/csv",
            key="codon_editor_csv_download",
        )
        render_cds_analysis_result(
            {
                "cds": edit["edited_cds"],
                "translated_amino_acids": edit["translated_amino_acids"],
                "analysis": last_result["analysis"],
            },
            title="替换后基础质检",
            key_prefix="codon_editor_result",
        )


def render_batch_tab(settings: dict[str, object]) -> None:
    uploaded = st.file_uploader("上传 AA FASTA 文件", type=["fa", "fasta", "faa", "txt"])
    raw_fasta = st.text_area("或粘贴 AA FASTA", value="", height=180)
    if st.button("批量预测", type="primary", key="batch_predict"):
        try:
            text = uploaded.read().decode("utf-8") if uploaded is not None else raw_fasta
            records = parse_fasta(text)
            results = []
            cds_records = []
            summary_rows = []
            for record in records:
                result = run_prediction(amino_acids=record.sequence, **settings)
                result["id"] = record.id
                result["description"] = record.description
                results.append(result)
                output_cds = result.get("postprocess", {}).get("optimized_cds", result["cds"])
                cds_records.append(FastaRecord(id=record.id, description="PichiaCLM optimized CDS", sequence=output_cds))
                analysis = result["analysis"]
                critical, warnings = quality_counts(analysis)
                summary_rows.append(
                    {
                        "ID": record.id,
                        "结论": quality_status(analysis),
                        "基础问题": critical,
                        "风险警告": warnings,
                        "AA 长度": len(result["amino_acids"]),
                        "CDS 长度": len(output_cds),
                        "GC%": analysis["gc_percent"],
                        "CAI 训练数据": analysis["cai"]["training"],
                        "CAI 公开表": analysis["cai"]["public"],
                        "翻译一致": analysis["translation_matches_input"],
                        "酶切位点": len(analysis["restriction_sites"]),
                        "Motif 命中": len(analysis["motif_hits"]),
                        "局部 GC 警告": len(analysis["local_gc_outliers"]),
                    }
                )
        except Exception as exc:
            st.error(str(exc))
            return

        st.subheader("批量结果")
        summary_rows.sort(key=lambda row: (row["基础问题"], row["风险警告"]), reverse=True)
        st.dataframe(summary_rows, use_container_width=True, hide_index=True)
        st.download_button("下载 CDS FASTA", data=format_fasta(cds_records), file_name="pichiaclm_batch_cds.fasta")
        st.download_button("下载结果表格 CSV", data=records_to_csv(summary_rows), file_name="pichiaclm_batch_report.csv")
        with st.expander("逐条详细结果"):
            for result in results:
                render_prediction_result(result, title=f"{result['id']} 详细结果", key_prefix=f"batch_{result['id']}")


def render_external_cds_tab(settings: dict[str, object]) -> None:
    st.subheader("二次优化 CDS 质检")
    st.caption("用于检查在外部软件中二次优化后的 CDS；这里不会重新调用模型预测。")

    single_cds = st.text_area("粘贴二次优化后的 CDS", value="", height=160, key="external_single_cds")
    expected_aa = st.text_area("可选：期望翻译得到的 AA", value="", height=120, key="external_expected_aa")
    if st.button("分析这条 CDS", type="primary", key="external_single_analyze"):
        try:
            result = run_cds_analysis(
                mode=settings["mode"],
                cds=single_cds,
                expected_amino_acids=expected_aa or None,
                api_url=settings["api_url"],
                motifs=settings["motifs"],
                custom_sites=settings["custom_sites"],
            )
        except Exception as exc:
            st.error(str(exc))
            return
        render_cds_analysis_result(result, key_prefix="external_single")

    st.divider()
    st.subheader("批量 CDS FASTA 质检")
    uploaded = st.file_uploader("上传 CDS FASTA 文件", type=["fa", "fasta", "fna", "txt"], key="external_cds_upload")
    raw_fasta = st.text_area("或粘贴 CDS FASTA", value="", height=180, key="external_cds_fasta")
    batch_expected_aa = st.text_area("可选：所有 CDS 共用的期望 AA", value="", height=100, key="external_batch_expected")

    if st.button("批量质检 CDS", type="primary", key="external_batch_analyze"):
        try:
            text = uploaded.read().decode("utf-8") if uploaded is not None else raw_fasta
            records = parse_fasta(text)
            results = []
            summary_rows = []
            for record in records:
                result = run_cds_analysis(
                    mode=settings["mode"],
                    cds=record.sequence,
                    expected_amino_acids=batch_expected_aa or None,
                    api_url=settings["api_url"],
                    motifs=settings["motifs"],
                    custom_sites=settings["custom_sites"],
                )
                result["id"] = record.id
                result["description"] = record.description
                results.append(result)
                analysis = result["analysis"]
                critical, warnings = quality_counts(analysis)
                summary_rows.append(
                    {
                        "ID": record.id,
                        "结论": quality_status(analysis),
                        "基础问题": critical,
                        "风险警告": warnings,
                        "CDS 长度": analysis["cds_length"],
                        "密码子数量": analysis["codon_count"],
                        "GC%": analysis["gc_percent"],
                        "CAI 训练数据": analysis["cai"]["training"],
                        "CAI 公开表": analysis["cai"]["public"],
                        "翻译一致": analysis["translation_matches_input"],
                        "酶切位点": len(analysis["restriction_sites"]),
                        "Motif 命中": len(analysis["motif_hits"]),
                        "局部 GC 警告": len(analysis["local_gc_outliers"]),
                        "非法碱基": ",".join(analysis["invalid_bases"]),
                    }
                )
        except Exception as exc:
            st.error(str(exc))
            return

        summary_rows.sort(key=lambda row: (row["基础问题"], row["风险警告"]), reverse=True)
        st.dataframe(summary_rows, use_container_width=True, hide_index=True)
        st.download_button("下载 CDS 质检 CSV", data=records_to_csv(summary_rows), file_name="pichiaclm_cds_qc.csv")
        st.download_button(
            "下载 CDS 质检 JSON",
            data=json_dumps({"records": results}),
            file_name="pichiaclm_cds_qc.json",
            mime="application/json",
        )
        with st.expander("逐条 CDS 质检详情"):
            for result in results:
                render_cds_analysis_result(
                    result,
                    title=f"{result['id']} 质检详情",
                    key_prefix=f"external_batch_{result['id']}",
                )


def render_fusion_tab(settings: dict[str, object]) -> None:
    signal = st.text_area("信号肽 AA", value="", height=100)
    mature = st.text_area("成熟蛋白 AA", value=DEFAULT_SEQUENCE, height=140)
    if st.button("比较整体优化与分段优化", type="primary", key="fusion_predict"):
        try:
            if settings["mode"] != DIRECT_MODE_LABEL:
                st.warning("信号肽拼接对比需要直接加载模型模式。")
                return
            predictor = load_predictor(settings["weights_path"], settings["device"])
            comparison = compare_signal_fusion(
                predictor,
                signal_peptide=signal,
                mature_protein=mature,
                allow_unknown=settings["allow_unknown"],
            )
            payload = asdict(comparison)
        except Exception as exc:
            st.error(str(exc))
            return

        st.subheader("信号肽拼接对比")
        st.metric("两种 CDS 是否完全一致", "是" if payload["cds_are_identical"] else "否")
        for label, key in [("整体优化", "whole_sequence"), ("分段优化", "segmented")]:
            item = payload[key]
            st.markdown(f"### {label}")
            render_prediction_result(
                {
                    **item["prediction"],
                    "analysis": item["analysis"],
                },
                title=f"{label}结果",
                key_prefix=f"fusion_{key}",
            )
            window = item["cleavage_window"]
            st.write(
                f"切割位点附近窗口: AA {window['amino_acid_start']}-{window['amino_acid_end']} / "
                f"CDS {window['cds_start']}-{window['cds_end']}"
            )
            st.code(window["amino_acids"], language="text")
            st.code(window["cds"], language="text")


def main() -> None:
    st.set_page_config(page_title="PichiaCLM", page_icon="DNA", layout="wide")
    st.title("PichiaCLM 氨基酸序列转 CDS")
    settings = sidebar_settings()
    tab_single, tab_candidates, tab_batch, tab_external, tab_editor, tab_fusion = st.tabs(
        ["单条预测", "多候选 CDS", "FASTA 批量", "二次优化 CDS 质检", "密码子编辑器", "信号肽拼接"]
    )
    with tab_single:
        render_single_tab(settings)
    with tab_candidates:
        render_candidates_tab(settings)
    with tab_batch:
        render_batch_tab(settings)
    with tab_external:
        render_external_cds_tab(settings)
    with tab_editor:
        render_codon_editor_tab(settings)
    with tab_fusion:
        render_fusion_tab(settings)


if __name__ == "__main__":
    main()

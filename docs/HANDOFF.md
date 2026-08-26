# PichiaCLM Handoff

## 当前切片

完成候选相似度硬门槛与 %MinMax 局部翻译速度谱两处修改：

- `core/candidates.py`：`select_low_similarity_subset` 新增显式 `max_bp_similarity_percent` / `max_codon_similarity_percent` 阈值，不满足时通过 `constraint_satisfied=False` 明确报告，不再静默返回近似解；阈值缺省时用 80/80，通过 `threshold_is_placeholder=True` 自报为非湿实验方确认的数字（机制见 ADR-0002；数值来源参考 Quan et al. 2011 同义变体库实测一致性，见 ADR-0004，取代 ADR-0002 中"纯占位无来源"的说法）。
- `core/analysis.py`：新增 `min_max_profile`（Clarke & Clark 2008 的 %MinMax 曲线）与 `compare_min_max_profiles`（两条曲线逐窗口比较）。
- `core/candidates.py`：新增 `temperature_sampling` 生成路径（与随机同义替换同级，ADR-0005）；`_rank_subset_by_min_max` 支持传入 `MinMaxHarmonizationTarget` 走 harmonization 排序（比对源物种曲线形状，ADR-0006），不传则沿用宿主最深负谷代理指标。
- `core/source_reference.py`：源物种密码子频率表与天然 CDS 的动态获取 + 本地缓存，支持手动输入，抓不到就报错不退回宿主数据（ADR-0006）。

以上都只在 `core` 层，未接入 `interfaces`（CLI/API/Streamlit）——研究人员目前只能通过 Python/notebook 直接调用。见 `docs/EXECUTION_PLAN.md`。

**当前最需要注意的一条**：相似度硬门槛的 bp 轴阈值（80%）在本流程下不可达——hLF 实测两两 bp 相似度 88.03%–93.80%，因为同义变体的密码子前两位被氨基酸锁死。结果是 `constraint_satisfied` 恒为 `False`，且永远因 bp 而失败。这需要用户决定怎么改（见 `docs/EXECUTION_PLAN.md` 待授权表），在此之前门槛的红色状态没有判别力。

**第二条**：`temperature_sampling` 的产出率取决于具体序列而非长度——hLF 前体（710 aa）跑满 10 条，更短的成熟链（691 aa）只出 4 条。代价是 CPU 上 364 秒 vs `kazusa_diverse` 的 3 秒，换来更好的多样性（最高密码子相似度 61.97% vs 72.39%）。实测表格见 `docs/EXECUTION_PLAN.md`。

另注意：`FakePredictor.predict_sample` 是**均匀采样**（`temperature` 在其中不起作用），与真实模型的 softmax 采样分布不同，不能用它估计这条路径的真实产出率。

## 下一步

优先解决上面那条已知限制（需要用户先决定 `temperature_sampling` 与 ADR-0001 硬门槛怎么共存）。其余待授权项：提供真实湿实验相似度阈值（替换 ADR-0004 的文献参考值）、把 `source_reference` 与 `MinMaxHarmonizationTarget` 的组装接到某一层、确定源基因 CDS 与设计的位点对齐口径、以及接口暴露。

## 必读材料

1. `docs/REQUIREMENTS.md`
2. `docs/ARCHITECTURE.md`
3. `docs/adr/ADR-0001-qualified-candidate-acceptance.md`
4. `docs/adr/ADR-0002-similarity-hard-threshold.md`
5. `docs/adr/ADR-0003-min-max-host-only-profile.md`
6. `docs/adr/ADR-0004-literature-informed-similarity-threshold.md`
7. `docs/adr/ADR-0005-temperature-sampling-strategy-and-min-max-ranking.md`
8. `docs/adr/ADR-0006-dynamic-source-reference-fetch.md`
9. `docs/EXECUTION_PLAN.md`

## 验证方式

运行 `python -m pytest -q tests/test_core_features.py tests/test_docs_governance.py`。若涉及接口、模型或部署，再增加与该切片对应的端到端验证。

## 硬约束

- PichiaCLM 输出是候选 CDS 和审查信息，不是表达产量、分泌效率或湿实验成功的预测。
- 不修改训练数据或模型权重，除非用户明确授权一个包含数据来源、评估门控和回滚方式的切片。
- 最终候选不得有关键问题、不得增加风险警告，也不得增加可避免的最低偏好密码子；CAI 不单独决定合格状态。
- 不改变远端可见性、不提交、不推送，除非用户明确要求。

```yaml
current_slice: candidate_similarity_gate_and_min_max_profile
slice_status: done
authorization_status: awaiting_user
verification_status: passed
```

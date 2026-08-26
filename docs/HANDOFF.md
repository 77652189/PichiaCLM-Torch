# PichiaCLM Handoff

## 当前切片

完成候选相似度硬门槛与 %MinMax 局部翻译速度谱两处修改，并已接入 CLI/API/Streamlit 三个接口：

- `core/candidates.py`：`select_low_similarity_subset` 新增显式 `max_codon_similarity_percent` 阈值，不满足时通过 `constraint_satisfied=False` 明确报告，不再静默返回近似解；缺省 80%，通过 `threshold_is_placeholder=True` 自报为非湿实验方确认的数字（机制见 ADR-0002；数值来源参考 Quan et al. 2011 同义变体库实测一致性，见 ADR-0004；只保留密码子轴见 ADR-0007）。
- `core/analysis.py`：新增 `min_max_profile`（Clarke & Clark 2008 的 %MinMax 曲线）与 `compare_min_max_profiles`（两条曲线逐窗口比较）。
- `core/candidates.py`：新增 `temperature_sampling` 生成路径（与随机同义替换同级，ADR-0005）；`_rank_subset_by_min_max` 支持传入 `MinMaxHarmonizationTarget` 走 harmonization 排序（比对源物种曲线形状，ADR-0006），不传则沿用宿主最深负谷代理指标。排序只重排 `selected_ranks`，不改变选中哪几条，也不与相似度门槛合成分数。
- `core/source_reference.py`：源物种密码子频率表与天然 CDS 的动态获取 + 本地缓存（缓存目录已 gitignore），支持手动输入，抓不到就报错不退回宿主数据（ADR-0006）。
- 接口层：`predict_candidates`、API 的 `PredictCandidatesRequest`、Streamlit 候选页都已打通 `strategy`、`max_codon_similarity_percent`、源物种 taxon 与天然 CDS。Streamlit 上门槛状态、排序依据、%MinMax 曲线与吻合度报告均可见。

**当前最需要注意的一条**：harmonization 目前只排序、不设计。候选生成阶段完全不知道目标曲线的存在，所以经常没有任何候选比基准设计更吻合源曲线——40 aa 测试序列实测最佳候选 9.11、基准 6.83。ADR-0009 的吻合度报告会把这种情况如实报出来（"没有候选优于基准"），不要把"排第一"读成"已完成 harmonization"。要真正拿到收益，需要让生成阶段朝目标曲线优化，这是待授权项。

**第二条**：`temperature_sampling` 的产出率取决于具体序列而非长度——hLF 前体（710 aa）跑满 10 条，更短的成熟链（691 aa）只出 4 条。决定因素是参考序列在 ADR-0001 门槛下的余量。代价是 CPU 上 364 秒 vs `kazusa_diverse` 的 3 秒，换来更好的多样性（最高密码子相似度 61.97% vs 72.39%）。实测表格见 `docs/EXECUTION_PLAN.md`。

**第三条**：相似度门槛的 80% 是文献参考值，不是湿实验方确认的上限，界面会自报占位。湿实验方给出数字后改一个参数即可，无需改代码逻辑。

另注意：`FakePredictor.predict_sample` 是**均匀采样**（`temperature` 在其中不起作用），与真实模型的 softmax 采样分布不同，不能用它估计这条路径的真实产出率。

## 下一步

无进行中的工作。所有剩余项都在 `docs/EXECUTION_PLAN.md` 的待授权表里，需要用户先做决定：让生成阶段朝目标曲线优化（优先，见上面第一条）、提供真实湿实验相似度阈值、决定 `temperature_sampling` 与 ADR-0001 硬门槛怎么共存、在 OPN（SPP1）上实测、把 `source_reference` 抓取结果接到调用方。

## 必读材料

1. `docs/REQUIREMENTS.md`
2. `docs/ARCHITECTURE.md`
3. `docs/adr/ADR-0001-qualified-candidate-acceptance.md`
4. `docs/adr/ADR-0002-similarity-hard-threshold.md`
5. `docs/adr/ADR-0003-min-max-host-only-profile.md`
6. `docs/adr/ADR-0004-literature-informed-similarity-threshold.md`
7. `docs/adr/ADR-0005-temperature-sampling-strategy-and-min-max-ranking.md`
8. `docs/adr/ADR-0006-dynamic-source-reference-fetch.md`
9. `docs/adr/ADR-0007-codon-axis-only-similarity-gate.md`
10. `docs/adr/ADR-0008-aligned-mature-peptide-source-cds.md`
11. `docs/adr/ADR-0009-harmonization-fit-against-baseline.md`
12. `docs/EXECUTION_PLAN.md`

## 验证方式

运行 `python -m pytest -q tests/`。其中 `tests/test_streamlit_candidates_tab.py` 用 Streamlit AppTest 驱动真实页面，覆盖接口层渲染（一次纯手工浏览器验证曾抓到只有整页渲染才暴露的 NameError，单元测试看不见）。

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

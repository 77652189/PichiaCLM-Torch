# PichiaCLM 执行计划

本文档是项目进度和授权范围的唯一权威。

## 当前状态

文档治理基线已建立：需求、架构、执行计划、handoff 与 ADR 索引各自承担单一职责。ADR-0001 记录现有最终候选判定的长期取舍。本轮由用户授权并完成一个具体科学切片（候选相似度硬门槛 + %MinMax 局部翻译速度谱，见 ADR-0002、ADR-0003），未授权新的模型、权重或部署行为。

## 已有能力摘要

- 已有核心预测、候选生成、CDS 分析、保守后处理和构建比较能力；
- 已有 CLI、FastAPI 与 Streamlit 三种接口；
- 最终候选相对基准的保守筛选与双 CAI 审查输出已存在；
- 候选子集筛选现在带显式相似度阈值门槛（未获湿实验数字前用文献参考值 80/80，自报为非湿实验方确认的占位值，见 ADR-0002、ADR-0004）。经代码核对：CLI/API/Streamlit 三个接口都通过 `core/predictor.py::PichiaCLMPredictor.predict_candidates` 这一个入口调用候选生成，该方法尚未新增 `max_bp_similarity_percent` / `max_codon_similarity_percent` 参数，三个接口也都没有对应的入参（CLI 无该 flag，API 的 `PredictCandidatesRequest` 无该字段，Streamlit 无对应输入框）——因此目前所有调用路径都只会用到 ADR-0002 的占位阈值，没有任何路径能传入真实数字。`constraint_satisfied` / `threshold_is_placeholder` 等新字段会随通用 `dict` 透传到 CLI 的 `--json` 输出和 API 响应体，但 CLI 的默认文本摘要和 Streamlit 界面都没有读取或展示这些字段（`interfaces/streamlit_app.py` 只读取 `selected_ranks`/`selected_size`/`requested_size`/`mean_codon_similarity_percent`/`max_codon_similarity_percent`）。
- `core/analysis.py` 新增 `min_max_profile` 与 `MinMaxWindow`（局部翻译速度曲线），以及 `compare_min_max_profiles` / `MinMaxProfileComparison`（两条曲线逐窗口比较，长度不一致直接报错，不做静默对齐）。
- `core/source_reference.py` 提供源物种参考数据的动态获取 + 本地缓存（人类密码子频率表、天然 CDS，支持手动输入），见 ADR-0006。
- `core/candidates.py` 的 `_rank_subset_by_min_max` 现在有两条排序口径：传入 `MinMaxHarmonizationTarget` 时按"候选在宿主频率下的曲线形状与源基因在源物种频率下的曲线的平均绝对差"排序（harmonization，越小越前）；不传时沿用原来的宿主最深负谷代理指标。两者都只重排 `selected_ranks`，不改变选中了哪几个候选，也不与相似度门槛合成分数（ADR-0002）。
- **已知限制（用真实权重实测，`Arch1-0404.weights.pt`，CPU）**：`temperature_sampling` 路径能否产出足够候选**取决于具体序列，不是简单的长度阈值**。实测（`num_candidates=10, seed=7, temperature=0.8`）：

  | 序列 | 长度 | 生成候选 | exhausted |
  | --- | ---: | ---: | --- |
  | 短控制序列 | 9 aa | 10 | 否 |
  | 白蛋白信号肽片段 | 60 aa | 10 | 否 |
  | 同上重复拼接 | 153 aa | 10 | 否 |
  | 另一种拼法 | 150 aa | 6 | **是** |
  | 同上延长 | 300 aa | 1（只有参考） | **是** |
  | **hLF 前体全长**（UniProt P02788） | **710 aa** | **10** | 否 |
  | **hLF 成熟链**（去 1–19 信号肽） | **691 aa** | **4** | **是** |

  hLF 实测确认这与长度无关：**更长的**前体（710 aa）跑满 10 条，**更短的**成熟链（691 aa）只出 4 条。决定因素是参考序列在 ADR-0001 门槛下的余量（前体 `avoidable_lowest=162`、成熟链 `151`），不是序列长度。

  代价差异很大：hLF 前体上 `kazusa_diverse` 用时 3 秒，`temperature_sampling` 用时 364 秒（CPU，约 120 倍）。但温度采样产出的多样性更好：推荐子集的最高密码子相似度 61.97% vs 72.39%。

  主导的拒绝原因是 ADR-0001 的既有硬门槛"可避免最低偏好密码子数不得高于参考"（150 aa 抽样 8 条全部因此被拒；300 aa 另有 3/8 因重复 k-mer 增加被拒）。参考序列在该门槛下留出的余量随序列不同而变化，所以 150 aa 会失败而 153 aa 不会。hLF（约 710 aa）单次运行在 CPU 上超过 2 分钟，本轮未跑完，按 300 aa 的趋势判断很可能退化为只有参考。

  兜底行为本身正确（明确报 `exhausted=True` 并给出 `note`，不静默）。

- **bp 相似度阈值 80% 在本流程下structurally 不可达（需用户定夺，ADR-0004 的遗留问题）**：同义变体必须保留氨基酸，密码子前两位基本被锁死，只有摇摆位可变，所以两条同义 CDS 的碱基相似度有一个远高于 80% 的下限。hLF 实测（`kazusa_diverse`，10 条候选，全部两两组合）：**bp 相似度 88.03%–93.80%，从未接近 80%**；同一组的密码子相似度是 69.3%–72.39%，轻松低于 80%。

  也就是说 `constraint_satisfied` 目前恒为 `False`，且失败原因永远来自 bp 那一侧，与候选质量无关。一个永远不可能通过的门槛会训练使用者忽略它——正是 `GUARDS.md` 里"家务型断言"要避免的形态。ADR-0004 把 Quan et al. 2011 的 79–82% 同时套到 bp 和 codon 两个轴上是错的：该数字用在 codon 轴上合理，用在 bp 轴上不可达。**需要用户决定**：bp 轴改用一个符合实际分布的上限、只保留 codon 轴作为硬门槛、还是维持现状并接受它恒红。在用户决定之前不擅自改数值。

- **测试替身与真实模型的采样分布不一致（需注意）**：`tests/test_core_features.py::FakePredictor.predict_sample` 用 `torch.ones(n) / temperature` 当权重，`multinomial` 归一化后是**均匀采样**，`temperature` 实际不起任何作用；真实模型走的是 `torch.softmax(logits / temperature)`，概率集中在偏好密码子上。用 FakePredictor 估计这条路径的产出率会得到明显偏悲观的结论（本轮曾据此得出错误结论，已用真实权重更正）。

## 待授权工作

| 工作 | 门控 | 现状 |
| --- | --- | --- |
| 用真实湿实验相似度阈值替换 ADR-0004 的文献参考值，并给 `predict_candidates` 及三个接口加对应入参 | 用户提供 bp/codon 相似度上限的具体数字（或湿实验方确认文献参考值可用） | 部分完成：数值来源已从"无来源占位"改为文献参考（ADR-0004），但仍非湿实验方确认，且三个接口都还没有入参 |
| 把候选池生成从随机同义位点替换换成 `predict_sample` 温度采样 | 用户确认这是要解决的问题，并给出验收标准 | 未授权 |
| 把 %MinMax 曲线接入至少一个接口展示 | 用户确认先接入哪个接口、以什么形式展示 | 未授权 |
| 提高 `temperature_sampling` 在长序列下的产出率（主因是 ADR-0001 的"可避免最低偏好密码子不劣于参考"门槛） | 用户确认怎么和 ADR-0001 共存：放宽门槛、改成按比例而非绝对数、允许按序列长度缩放、还是接受现状 | 未授权。该路径在部分序列上可用、部分不可用，取决于参考序列在该门槛下的余量 |
| 修正 bp 相似度阈值（当前 80% 不可达，`constraint_satisfied` 恒为 False） | 用户决定：bp 轴换成符合实测分布的上限 / 只用 codon 轴做硬门槛 / 维持现状 | **未授权，阻塞硬门槛的实际可用性** |
| ~~在 hLF 真实序列上实测~~ | — | **已完成**：UniProt P02788，前体 710 aa 与成熟链 691 aa 均已测，数据见上 |
| 在 OPN（SPP1）真实序列上实测 | 用户确认使用哪个转录变体 | 未验证 |
| 把 `source_reference` 的抓取结果接到实际调用方（构造 `MinMaxHarmonizationTarget` 的那一步现在要调用方自己做） | 用户确认在哪一层组装：`predict_candidates` 参数、CLI/API 字段，还是留给 notebook | 未授权 |
| 校验源基因 CDS 与候选设计的位点对齐（`compare_min_max_profiles` 要求两条曲线窗口数相同；人类天然 CDS 含信号肽/终止密码子，未必与设计逐位对应） | 用户确认对齐口径：截掉信号肽、按成熟肽对齐、还是要求研发组同事输入已对齐的 CDS | 未授权 |
| 改变模型、权重或训练数据 | 用户确认数据来源、评估指标及回滚方式 | 未授权 |
| 改变 LAN/API 部署行为 | 用户确认目标环境与可见性范围 | 未授权 |

## 明确不做

- 不因文档更新而修改现有模型、推理逻辑、接口行为或部署配置；
- 不将历史 UI 表格、候选分数或本地服务状态当作当前验收证据；它们必须在相应切片中重新验证。

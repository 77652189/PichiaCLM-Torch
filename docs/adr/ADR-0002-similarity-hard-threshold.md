# ADR-0002：候选相似度筛选改为显式阈值门槛，不做静默最小化

**状态：** accepted<br>
**日期：** 2026-08-25

## 背景

湿实验负责人要求：送去做实验的候选序列之间不能相似度过高。此前 `select_low_similarity_subset` 只是一个最小化器：在给定子集大小下挑两两相似度总和最低的一组，没有阈值参数。如果所有候选都过于相似，它仍然返回"最不坏"的一组，且不报告约束是否满足——不满足的情况要等送样被拒才暴露。

此前制造多样性的做法是对参考 CDS 随机改动约 10%–20% 的同义位点（见 `core/candidates.py` 的 `min_difference_percent` / `max_difference_percent`）。这个做法本身也有问题：随机同义替换是靠破坏结构换取多样性——每一次随机替换，抹平天然翻译暂停位点的概率与保留它的概率相当，而局部翻译速度分布正是共翻译折叠所依赖的东西。多样性更应该来自 `core/predictor.py` 的 `predict_sample` 温度采样。但把候选池生成机制从随机同义替换换成温度采样是更大的改动，需要单独授权的切片，不在本次修改范围内；本 ADR 只处理"如何判定一组候选是否够不相似"。

## 决策

给 `select_low_similarity_subset` 和 `CandidateGenerationOptions` 加两个显式的最大相似度参数：`max_bp_similarity_percent`、`max_codon_similarity_percent`（`CandidatePairSimilarity` 里已经分别有 `bp_similarity_percent` 和 `codon_similarity_percent` 两种同一性）。

返回结果 `CandidateSubsetSelection` 新增：

- `bp_similarity_threshold_percent` / `codon_similarity_threshold_percent`：本次实际使用的阈值；
- `threshold_is_placeholder`：阈值是否未由调用方显式提供；
- `constraint_satisfied`：挑出的子集（仍然是两两相似度最低的那组）是否满足阈值。

不满足阈值时，函数仍然返回这组"最不坏"的候选（而不是 `None` 或抛异常），但 `constraint_satisfied=False`，调用方可以直接从已有的 `max_bp_similarity_percent` / `max_codon_similarity_percent` 与对应阈值字段算出差多少，不需要额外的"差值"字段。

阈值数值来自湿实验方，不是能从代码推导的。在拿到具体数字之前，两个阈值都使用占位值 `PLACEHOLDER_MAX_BP_SIMILARITY_PERCENT` / `PLACEHOLDER_MAX_CODON_SIMILARITY_PERCENT`（当前均为 80.0），并通过 `threshold_is_placeholder=True` 自报为占位值，不得被解读为生物学判断。调用方一旦获得湿实验的真实阈值，通过 `CandidateGenerationOptions.max_bp_similarity_percent` / `max_codon_similarity_percent` 显式传入即可关闭占位标记。

## 后果

- `select_low_similarity_subset` 默认就带门槛（占位值），而不是"不设阈值"作为常态默认——避免这一检查被当作可选项而遗忘配置，退回到静默返回近似解。
- CLI / API / Streamlit 在未显式配置真实阈值前看到的 `constraint_satisfied` 都是基于占位值的诊断，不能展示为已获批准的相似度上限。
- 候选池生成仍依赖随机同义位点替换制造多样性，而不是温度采样；这是背景里指出的问题根源，但换成温度采样需要单独授权的切片，记入 `docs/HANDOFF.md` 的下一步。

## 备选方案

- **不加阈值，只做最小化（维持现状）**：拒绝，因为约束是否满足会一直推迟到送样被拒才暴露，违反"不得用空表、默认值或伪造记录静默继续"的既有原则。
- **把相似度和 %MinMax 吻合度（见 ADR-0003）合成一个综合分再排序**：拒绝，因为合成分会让硬约束被其他项的高分补偿掉，可能交付一组根本不满足湿实验要求的序列。两两相似度是约束（可以失败），%MinMax 吻合度是偏好（只在通过门槛的候选之间排序），不是同一回事，不能用同一个数表达。
- **继续依赖随机改变约 10% 同义位点来获得多样性，不引入独立的相似度判定**：拒绝，因为这是靠破坏结构换取多样性，且不检验结果是否真的低相似度。

## 取代关系

关于阈值数值来源的说明（"阈值数值来自湿实验方，不是能从代码推导的"一句）已被 [ADR-0004](ADR-0004-literature-informed-similarity-threshold.md) 取代：数值不变，但改为参考 Quan et al. 2011 的文献先例，而非纯粹无来源的占位值。本 ADR 关于硬门槛机制本身（`constraint_satisfied`、不满足时不静默返回近似解、`CandidateSubsetSelection` 新增字段等）的决策未被取代，继续有效。

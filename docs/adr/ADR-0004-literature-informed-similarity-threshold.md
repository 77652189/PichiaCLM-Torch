# ADR-0004：相似度阈值改为文献参考值，取代 ADR-0002 中"纯占位、无来源"的说法

**状态：** accepted（bp 轴部分已被 [ADR-0007](ADR-0007-codon-axis-only-similarity-gate.md) 取代）<br>
**日期：** 2026-08-25<br>
**取代：** [ADR-0002](ADR-0002-similarity-hard-threshold.md) 中关于阈值数值来源的说明（不影响 ADR-0002 关于硬门槛机制本身的决策）

## 背景

ADR-0002 把 `max_bp_similarity_percent` / `max_codon_similarity_percent` 都设为占位值 80.0，并声明"阈值数值来自湿实验方，不是能从代码推导的……不得被解读为生物学判断"——也就是纯粹任意的占位数字，没有任何来源。

用户要求：在没有湿实验方给出具体数字之前，先从湿实验的需求做合理推断，或联网查文献，而不是继续用一个毫无来源的数字。

查到的最贴近先例是 Quan et al. 2011（*Nature Biotechnology*，"Parallel on-chip gene synthesis and application to optimization of protein expression"）：他们为同一蛋白（scFv、polymerase）构建纯同义密码子变体库（scFv 21 个变体、polymerase 24 个变体），通过 Monte Carlo 采样让变体尽量分散，目的是并行合成后做表达筛选——这与本项目"生成多个候选 CDS 送湿实验筛选"是同一类问题。该库最终库内两两 DNA 序列一致性的**平均值**为 79%（scFv）、82%（polymerase）。

这个数字不能直接照抄，因为统计口径不同：

- Quan et al. 报告的是"一个大库（20+ 变体）里所有两两组合的**平均**一致性"；
- 本仓库的门槛判的是"选中子集（通常 3–5 个）里最差一对的**最大**一致性"（`max_bp_similarity_percent` / `max_codon_similarity_percent`）。

对同一组序列，两两相似度的最大值天然大于或等于平均值。所以如果把"80%"直接当作 max 上限，标准并不比 Quan et al. 实际做到的更松——即便统计口径不同，用同一个数字作为 max 上限至少不会比文献先例更宽松。

另外还查到 JigsawSeq（Nauman et al. 2015，*Nature Communications*）等方法：用物理连接在同一分子上的分子条形码区分变体身份，CDS 主体本身不需要强制拉低相似度。如果湿实验要求低相似度的真实动机是"区分测序/PCR 结果里哪个候选是哪个"，条形码可能是更标准的解法；但这是与湿实验协议相关的选择，属于用户与湿实验方之间的决定，本 ADR 不代为决定，只记录这是一个存在的替代路径。

## 决策

阈值数值维持 80.0（bp）/ 80.0（codon）不变，但来源说明改变：

- 不再声称是"纯粹任意、无来源"的占位值；改为"参考 Quan et al. 2011 对同义变体库的实测两两一致性（79%–82% 均值），取整近似为 80%，作为 max 上限使用"。
- `threshold_is_placeholder=True` 的语义保持不变：这仍然不是湿实验方直接确认的数字，调用方看到这个标记时仍不能把它当作已获批准的相似度上限。
- `core/candidates.py` 中 `PLACEHOLDER_MAX_BP_SIMILARITY_PERCENT` / `PLACEHOLDER_MAX_CODON_SIMILARITY_PERCENT` 的注释同步更新为指向本 ADR，说明数值来源与统计口径差异，而不是继续说"没有来源"。

## 后果

- 阈值数值本身没有变化（仍是 80.0 / 80.0），所以现有代码行为和测试不受影响；这是一次"说明来源"的变更，不是"改变行为"的变更。
- 一旦湿实验方给出真实数字，替换方式不变：通过 `CandidateGenerationOptions.max_bp_similarity_percent` / `max_codon_similarity_percent` 显式传入，`threshold_is_placeholder` 会随之变为 `False`。
- `docs/EXECUTION_PLAN.md` 的待授权工作里，"用真实湿实验相似度阈值替换占位值"这一项状态从"未授权"改为"临时用文献参考值替代，仍待湿实验方确认或替换"。
- 条形码类替代方案未被采纳，也未被排除；记入 `docs/HANDOFF.md`，供用户后续与湿实验方确认协议本身的需求时参考。

## 备选方案

- **继续维持 ADR-0002 的"纯占位、无来源"说法，不做任何调整**：拒绝，因为用户明确要求先做合理推断，而不是放着一个毫无依据的数字。
- **直接把 Quan et al. 的 79%/82% 平均值原样当作 max 上限，不说明统计口径差异**：拒绝，因为均值和最大值是不同的统计量，不说明清楚会让读者误以为这是同口径的直接引用，构成"未经真实数据验证的科学结论"式的误导。
- **改用条形码方案替代整条 CDS 去相似化**：拒绝（本 ADR 范围内不采纳），因为这是协议层面的改动，需要用户和湿实验方确认这是否解决的是同一个问题；本次只记录该路径存在，不代为决定。

## 取代关系

取代 [ADR-0002](ADR-0002-similarity-hard-threshold.md) 中"阈值数值来自湿实验方，不是能从代码推导的"这一句关于来源的说法；ADR-0002 关于硬门槛机制（`constraint_satisfied`、不满足时不静默返回近似解等）的决策继续有效，未被取代。

**本 ADR 中"两个阈值都取 80.0"里关于 bp 轴的部分已被 [ADR-0007](ADR-0007-codon-axis-only-similarity-gate.md) 取代**：hLF 实测显示同义变体的 bp 相似度下限远高于 80%，该轴的门槛不可达，现已移除，bp 相似度降级为展示信息。本 ADR 关于 codon 轴取值 80.0、来源为 Quan et al. 2011、以及仍属占位值的定性未被取代，继续有效。

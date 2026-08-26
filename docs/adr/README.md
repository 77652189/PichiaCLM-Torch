# 架构决策记录索引

| 编号 | 标题 | 状态 | 摘要 |
| --- | --- | --- | --- |
| [ADR-0001](ADR-0001-qualified-candidate-acceptance.md) | 最终候选采用相对基准的保守筛选 | accepted | 以质量风险不劣化为合格门槛；CAI 提供双参考审查，不作为单独阈值。 |
| [ADR-0002](ADR-0002-similarity-hard-threshold.md) | 候选相似度筛选改为显式阈值门槛 | accepted（阈值来源说明已被 ADR-0004 取代） | 相似度门槛可显式设置，不满足时明确报告而非静默返回近似解；阈值缺省时使用自报的占位值。 |
| [ADR-0003](ADR-0003-min-max-host-only-profile.md) | %MinMax 局部翻译速度谱先算宿主曲线 | accepted（"源物种因缺数据推迟"的框定已被 ADR-0005 取代） | 当前只计算相对宿主频率的曲线。 |
| [ADR-0004](ADR-0004-literature-informed-similarity-threshold.md) | 相似度阈值改为文献参考值 | accepted（bp 轴部分已被 ADR-0007 取代） | 阈值来源改为参考 Quan et al. 2011 同义变体库实测一致性，而非无来源占位值；仍非湿实验方确认的数字。 |
| [ADR-0005](ADR-0005-temperature-sampling-strategy-and-min-max-ranking.md) | 新增温度采样策略，用 %MinMax 排序其结果 | accepted（"源物种对比不适用"一条已被 ADR-0006 取代） | 与随机同义替换同级的第二条生成路径，不套用 10%–20% 差异区间；%MinMax 只排序通过相似度门槛的候选，不合成分数。 |
| [ADR-0006](ADR-0006-dynamic-source-reference-fetch.md) | 源物种参考数据动态获取 + 本地缓存 | accepted | hLF/OPN 是人源异源表达，源物种对比适用；人类密码子频率表与天然 CDS 都按需抓取并缓存，支持手动输入，抓不到就报错不退回宿主数据。 |
| [ADR-0007](ADR-0007-codon-axis-only-similarity-gate.md) | 硬门槛只用密码子轴 | accepted | hLF 实测 bp 相似度下限 88%，80% 阈值不可达导致门槛恒红；bp 轴移除，降级为展示信息，codon 轴照常。 |

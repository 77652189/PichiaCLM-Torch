# ADR-0006：源物种参考数据用动态获取 + 本地缓存，不硬编码进代码

**状态：** accepted<br>
**日期：** 2026-08-26<br>
**取代：** [ADR-0005](ADR-0005-temperature-sampling-strategy-and-min-max-ranking.md) 中"源物种对比对本项目不适用"的结论（不影响 ADR-0005 关于温度采样生成路径、以及 %MinMax 只排序不合成分数的决策）

## 背景

ADR-0003 记录过：本仓库现在的 `min_max_profile` 只能算相对宿主（毕赤酵母）的曲线，因为没有源物种密码子使用表。ADR-0005 随后把这件事进一步框定为"源物种对比对本项目不适用"，依据是 PichiaCLM 的优化对象、训练数据、湿实验宿主三者都是毕赤酵母。

用户现在澄清了一个当时没说明的事实：一兮生物实习项目里要表达的 hLF、OPN 是**人源蛋白**（源物种 Homo sapiens，taxon 9606），在毕赤酵母里做异源表达。这正是文献里 codon harmonization（Wright et al. 2022，CHARMING）针对的场景——宿主是毕赤酵母不假，但蛋白本身有一个不同的源物种，所以两条曲线都存在。ADR-0005 的"不适用"结论建立在"三者同物种"这个前提上，该前提对这两个目标蛋白不成立，因此那条结论需要更正。

在此基础上进一步确认：只有人类密码子使用频率表还不够。%MinMax 要比较的是"某个位置历史上实际用的密码子，相对源物种整体频率，是常用还是罕见"——这是一个**密码子层面**的事实，氨基酸序列是简并的，同一个氨基酸在源物种里可能被编码成任意一个同义密码子，氨基酸序列本身不携带"当年具体选了哪一个"这个信息。所以要重建源物种的 %MinMax 目标曲线，除了频率表，还需要 hLF（基因 LTF）、OPN（基因 SPP1）**真实的人类天然 CDS 核苷酸序列**本身。

用户明确要求：这两类参考数据（人类密码子频率表、天然 CDS）都不应该硬编码进代码，应该做成"动态获取 + 本地缓存"，并且天然 CDS 还要支持研究人员手动输入——因为 SPP1 本身有多个转录变体，自动抓到的默认版本不一定是研究人员实际在用的那个构建。

## 决策

新增 `core/source_reference.py`，提供两个函数，两者都遵循同一套解析顺序：**本地缓存 → 外部数据源 → 失败就报错**，不静默退回宿主数据或空表：

1. `load_source_organism_codon_fractions(taxon_id, *, cache_dir=None, fetch=..., timeout=10.0) -> (fractions, total_codon_count)`
   - 数据源：Kazusa Codon Usage Database（`showcodon.cgi?species={taxon_id}&aa=1&style=N`）——跟本模块里已经手工抄录的宿主表 `PUBLIC_PICHIA_PASTORIS_FRACTIONS` 同一个数据库、同一种查询方式，口径一致。
   - 不加 `python_codon_tables` 这类第三方库依赖：Kazusa 的返回是纯文本表，几十行代码能自己解析，没必要为此引入新依赖。
2. `load_native_source_cds(*, accession=None, manual_cds=None, cache_dir=None, fetch=..., timeout=10.0) -> str`
   - `manual_cds` 优先于一切：给了就直接用（本地校验+标准化），完全不碰网络、也不缓存——这是研究人员用来确认"就是这个构建"的直接通道。
   - 否则按 `accession` 走 NCBI RefSeq 的 `efetch`（`rettype=fasta_cds_na`）。

两个函数都接受可注入的 `fetch` 参数（默认是基于 `urllib.request` 的实现），测试和未来任何调用方都可以换成假的抓取函数，不需要真的连网络。

**缓存**：本地 JSON 文件，默认存在 `Model_PichiaCLM/Training/ExternalReferenceCache/`，按 `taxon_id` 或 `accession` 命名。这个目录**不进 git**（已加进 `.gitignore`）：这是"抓来的事实"，不是像 `Training/AllData` 那样经过校验、固化的训练资产，换一台机器重新抓一次没有成本。缓存没有 TTL 自动过期；要刷新就删缓存文件。

**失败时的兜底**：网络失败且本地无缓存 → 抛 `RuntimeError`，说明抓不到、且没有回退到宿主数据或空表。`manual_cds` 格式不对（含非 ATGCU 字符）→ 抛 `ValueError`。这跟 ADR-0002/ADR-0003 一贯的"不静默"原则一致：源物种数据抓不到，就是这个功能现在不可用，不是安静地退化成别的东西。

## 后果

- `min_max_profile` 本身不用改：它已经是"给什么频率表就用什么表"的通用签名（ADR-0003 时就是这么设计的），`load_source_organism_codon_fractions` 的返回值可以直接喂给它。
- 真正的 harmonization 排序（拿源物种曲线当目标，给候选排序）**还没接线**：这次只解决"怎么把源物种的两类参考数据弄到手"，不涉及"怎么用这两类数据去给候选排序"——那需要另外定义"曲线吻合度"怎么算，留在 `EXECUTION_PLAN.md` 待授权工作里。
- 离线环境下，任何依赖这两个函数的功能都会直接报错而不是退化，调用方需要能处理 `RuntimeError`（比如提示用户改用 `manual_cds` 或先联网跑一次把缓存建好）。
- 新增了对 `kazusa.or.jp` 和 `eutils.ncbi.nlm.nih.gov` 两个外部服务的网络依赖（仅在缓存未命中时触发）；`core` 包顶部的注释"intentionally has no dependency on FastAPI, Streamlit, or CLI code"指的是不依赖这些框架，不是不能有任何 I/O，这里没有违反——但这是 `core` 里第一次出现网络调用，值得记录在案。

## 备选方案

- **把人类密码子频率表和 hLF/OPN 的 CDS 直接抄进代码**（像现在的 `PUBLIC_PICHIA_PASTORIS_FRACTIONS` 一样）：拒绝，用户明确要求不要硬编码；而且天然 CDS 会随着"研究人员实际在用哪个构建"而变，硬编码一份没法覆盖 SPP1 多个转录变体的情况。
- **加 `python_codon_tables` 第三方库来做频率表抓取**：拒绝，Kazusa 的纯文本格式足够简单，没必要为此新增一个运行时依赖；且该库不内置人类表，仍然要走它自己的网络抓取路径，没有省掉网络依赖，只是换了一层。
- **只支持自动抓取天然 CDS，不支持手动输入**：拒绝，用户明确要求要支持手动输入；且 SPP1 有多个转录变体，自动抓的默认 accession 不一定是研究人员实际在用的序列，手动输入是必要的精度保障，不是可选的便利功能。
- **抓取失败时退回宿主频率表凑合用**：拒绝，这正是 ADR-0003 已经明确否定过的"用宿主数据冒充源物种数据"，会让 harmonization 目标曲线看起来存在但其实没有意义。

## 取代关系

取代 [ADR-0005](ADR-0005-temperature-sampling-strategy-and-min-max-ranking.md) 的"源物种对比对本项目不适用、不是数据缺口"这一条结论：对 hLF/OPN 这类人源异源表达目标，源物种对比适用，需要的源物种数据由本 ADR 的动态获取 + 缓存机制提供。ADR-0005 的其余决策（温度采样生成路径、%MinMax 只排序不过滤不合成分数、接口暴露范围）未被取代，继续有效。

本 ADR 如被后续决策取代，应新建 ADR 并在本文件及索引中标明替代关系。

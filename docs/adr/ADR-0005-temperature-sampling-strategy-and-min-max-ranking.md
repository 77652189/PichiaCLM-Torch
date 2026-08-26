# ADR-0005：新增温度采样候选生成路径，并用 %MinMax 排序其结果；源物种对比对本项目不适用

**状态：** accepted（"源物种对比不适用"这一条结论已被 [ADR-0006](ADR-0006-dynamic-source-reference-fetch.md) 取代）<br>
**日期：** 2026-08-25<br>
**取代：** [ADR-0003](ADR-0003-min-max-host-only-profile.md) 中"源物种对比因缺数据被推迟"的框定（不影响 ADR-0003 关于 `min_max_profile` 只算宿主曲线这一实现决策本身）

## 背景

现有候选生成只有一条路径：`core/candidates.py` 的 `_generate_candidate_pool` 对参考 CDS 随机改动 10%–20% 的同义位点。用户已经指出这个路径本身有问题——随机同义替换是"靠破坏结构换多样性"，多样性应该来自 `core/predictor.py::predict_sample` 的温度采样。

用户要求：加一条同级别的新路径（温度采样），让研究人员自己选，并把 %MinMax 接进去。

ADR-0003 之前把"源物种曲线对比"框定为"更有意义、但因为仓库里没有源物种密码子频率表而推迟"。用户澄清：**本项目需要的就是毕赤酵母（宿主）视角，训练数据本身也是毕赤酵母的**——也就是说，源物种对比不是"数据还没到位"，而是这个项目根本不需要它：PichiaCLM 要优化的对象、训练模型用的数据、湿实验要验证的宿主，三者都是同一个物种。ADR-0003 的框定需要更正。

## 决策

### 1. 源物种对比：不适用，不是推迟

`min_max_profile` 继续只算相对某个频率表的曲线（不变）；但不再把"缺源物种数据"当作一个待补的缺口。对 PichiaCLM 而言，宿主曲线（用训练数据频率表算，因为训练数据就是毕赤酵母）就是唯一有意义的比较对象。`docs/EXECUTION_PLAN.md` 里"接入源物种密码子使用频率表"这一条待授权工作按此撤销。

### 2. 新增 `temperature_sampling` 生成策略

`CandidateGenerationOptions.strategy` 新增合法值 `"temperature_sampling"`（原有 `"kazusa_diverse"` 保持默认，行为不变）。`_generate_candidate_pool` 按 `strategy` 分派：

- `"kazusa_diverse"`：走原有逻辑，不变。
- `"temperature_sampling"`：反复调用 `predictor.predict_sample(reference.amino_acids, temperature=options.temperature, generator=<按 options.seed 播种的 torch.Generator>)`，每次采样得到一条完整候选 CDS，和参考 CDS 逐密码子比较得到 `changes`，复用现有的 `CandidateDraft`/去重签名/质量过滤逻辑。

两条路径共享、不重复实现的部分：翻译一致性检查、GC/homopolymer/重复/motif/限制性位点风险不劣于参考、可避免最低偏好密码子数不劣于参考（ADR-0001 的既有硬约束，与生成策略无关）、之后的 `_pop_most_diverse_draft` 贪心多样性挑选、`select_low_similarity_subset` 相似度硬门槛（ADR-0002/0004）。

**`temperature_sampling` 路径不套用 `min_difference_percent`/`max_difference_percent`（10%–20%）区间过滤**：这个区间正是"随机改约 10% 同义位点"那条被质疑的机制的直接产物；温度采样的多样性应该由采样过程本身决定，继续用同一个区间卡它，等于换了个生成器却留着旧生成器的强制约束，没有真正解决问题。

### 3. %MinMax 只用来给通过相似度门槛的候选排序，不合成分数，不做过滤

只在 `strategy == "temperature_sampling"` 时生效：`recommended_subset.selected_ranks` 里的候选，按各自 CDS 相对训练数据频率表算出的 %MinMax 曲线中**最深的负值谷**（`min(windows.percent)`，忽略 `percent=None` 的窗口）重新排序——谷越浅（越不负）排越前。这只是排序，不影响 `select_low_similarity_subset` 已经选出的是哪几个候选，也不产生一个把相似度和 %MinMax 揉在一起的综合分（沿用 ADR-0002 的原则：约束和偏好不能用同一个数表达）。

无法算出谷值的候选（比如 CDS 长度不足一个窗口）排在最后，不当作"最差"处理，因为那是"没法比较"不是"比较后更差"。

### 4. 接口暴露范围

这一轮只把 `strategy` 参数打通到 `core/predictor.py::PichiaCLMPredictor.predict_candidates`（`CandidateGenerationOptions` 本来就有这个字段，`predict_candidates` 之前没有转发它）。CLI 的 `--strategy` flag、API 的 `PredictCandidatesRequest.strategy` 字段、Streamlit 的策略选择控件都没有加——这跟 ADR-0002/0003 的先例一致：核心能力先落地、接口暴露留给用户后续单独授权。也就是说"研究人员自行选择"目前是指"能通过 Python/notebook 直接调用 `predict_candidates(strategy=...)`"，还不是"能在网页界面上点选"。

## 后果

- 新策略生成的候选仍然受 ADR-0001 的合格判定约束（不劣于参考的质量/风险/可避免低偏好密码子数），不会因为换了生成方式就放宽这些硬约束。
- 温度采样路径下，如果模型在给定温度下总是采样出和参考几乎一样的结果（`changes` 为空），会被当作"没有产生新候选"跳过，可能导致池子填不满、最终 `exhausted=True` 并在 `note` 里报告——这是既有的兜底机制，不是本次新增的行为。
- %MinMax 排序依赖 `load_training_codon_reference()`，跟现有 CAI 的"训练数据参考"用的是同一份数据源，没有引入新的数据依赖。
- CLI/API/Streamlit 还是看不到 `strategy` 选项，也看不到 %MinMax 排序的效果，除非用户单独授权接口暴露这一步。

## 备选方案

- **继续对温度采样路径套用 10%–20% 差异区间过滤**：拒绝，理由见上——这会保留被质疑的那个约束机制，没有真正换掉生成方式背后的假设。
- **把相似度和 %MinMax 谷值合成一个分数排序**：拒绝，沿用 ADR-0002 已经否掉的理由：约束可以失败、偏好只用来排序，合成一个数会让约束被偏好项的高分掩盖。
- **等有源物种数据再做 %MinMax，本次不接入排序**：拒绝，用户已经明确"我们需要的就是毕赤酵母的"——不存在"等数据"这件事，宿主曲线就是最终要的东西，没有理由再往后拖。
- **直接把 strategy 也接进 CLI/API/Streamlit**：这次没做，不是因为不该做，是因为没必要在同一轮里把核心能力和三个接口的改动都做完；接口暴露范围更小、更容易独立验证，留给用户明确要不要现在做。

## 取代关系

取代 [ADR-0003](ADR-0003-min-max-host-only-profile.md) 中"源物种对比是更有意义的比较、只是因为缺数据而推迟"的框定。ADR-0003 关于 `min_max_profile` 本身只算宿主频率曲线的实现决策未被取代，继续有效。

**本 ADR 的第 1 节"源物种对比：不适用，不是推迟"已被 [ADR-0006](ADR-0006-dynamic-source-reference-fetch.md) 取代**：当时的依据是"优化对象、训练数据、湿实验宿主三者同为毕赤酵母"，但 hLF/OPN 是人源蛋白在毕赤酵母里的异源表达，该前提不成立。本 ADR 其余各节（温度采样生成路径、%MinMax 只排序不过滤不合成分数、接口暴露范围）未被取代，继续有效。

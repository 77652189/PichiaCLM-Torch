# PichiaCLM 密码子设计与可构建候选

[English](README.md) · [中文](README.zh.md)

> 把一条蛋白序列变成若干条**彼此互不相同**的同义 CDS 候选，
> 每条都过序列层面的构建风险筛查，且在关键指标上都不许比基准更差。

"彼此互不相同"是整件事的要害。一个返回五条、实际上是同一条序列的生成器，
返回的是**一条**候选——湿实验没法拿它当五个构建体去跑。

---

## 贡献边界

本仓库 fork 自
[owen-min/PichiaCLM-Torch](https://github.com/owen-min/PichiaCLM-Torch)。

**模型的 Keras → PyTorch 移植不是我做的。** 那是上游基线，提交 `b6fea1b`。
我的贡献是它之后的全部：序列安全校验、带约束的多候选生成、两两互异的子集选择、
密码子编辑器，以及三个共用同一核心的入口。`git log` 里这条分界线是直接可查的。

## 模型

不是 Transformer——是**多任务 GRU Seq2Seq + 缩放点积注意力**
（[`core/model.py`](Model_PichiaCLM/core/model.py)）：

- 在氨基酸序列上做双向 GRU 编码
- 一个 GRU 解码器输出密码子，另一个 GRU 解码器把氨基酸序列重建出来作为辅助任务——
  这个重建头的作用是把密码子头**拴在**它本该编码的那条蛋白上
- 解码状态与编码输出之间做缩放点积注意力（`dot_product_attention`）

权重随仓库分发（约 37 MB），推理在 CPU 上跑。不需要 GPU，不需要下载步骤。

## 方法：一条候选怎样才算合格

采样同义密码子很容易。难的是让候选**实验室真能构建**。

**① 带约束的解码。** 解码时做掩码，只允许输出与目标残基同义的密码子。
于是"翻译回去还是同一条蛋白"是结构上保证的，
而不是事后再检查一遍、然后祈祷它通过。

**② 序列安全校验**（[`core/analysis.py`](Model_PichiaCLM/core/analysis.py)）。
每条候选都要过一遍真正会让合成与克隆失败的那些模式——
而不是那些在指标上好看的模式：

| 检查项 | 为什么要有 |
| --- | --- |
| `LocalGCWindow` | 全局 GC 正常时，局部 GC 极端照样让合成失败 |
| `HomopolymerRun` | 单碱基长串导致聚合酶打滑 |
| `TandemRepeat` · `RepeatedKmer` | 重复序列导致组装错位 |
| `MotifHit` | 非预期的限制性酶切位点与调控元件 |
| `CAIComparison` | 密码子适应度，对着**两套**参照系各算一遍 |

**③ 两套参照系，而且都不当闸门。** CAI 同时对训练集密码子频率和公开 Kazusa 频率各算一次。
两者会不一致——而这个不一致本身就是信息，
所以两个值并排列在候选旁边，谁都不是通过/不通过的阈值
（[ADR-0001](docs/adr/ADR-0001-qualified-candidate-acceptance.md)）。
只设一个 CAI 阈值会简单得多，代价是把"这个判决出自哪套参照系"藏起来了。

**④ 合格是相对基准的，不是绝对的。** 一条新候选被接受，必须同时满足：
翻译回同一条蛋白、没有关键质量问题、风险警告数**不高于**基准、
可避免的最低偏好密码子数**不高于**基准。
CAI 更高的候选照样可能被拒——
用增加的风险去换一个更好看的分数，正是这条规则要挡住的交易。

**⑤ 多样性是选出来的，不是指望出来的**
（[`core/candidates.py`](Model_PichiaCLM/core/candidates.py)）。
`PairwiseDiversity` 度量候选之间在碱基与密码子层面差多远，
`CandidateSubsetSelection` 挑出相互差异最大的那个子集。
设计空间耗尽时，`CandidateSet.exhausted` 置位，**返回的候选就是少于请求数**——
唯一绝对不能做的，是拿近似重复的序列把数量凑满。

最后这条最值得在评审里守住：请求 5 条返回 3 条，是一个看得见、可以补救的失望；
返回 5 条、其实只有 3 条，会在任何人察觉之前先烧掉几周湿实验。

## 快速开始

```powershell
python -m streamlit run Model_PichiaCLM/interfaces/streamlit_app.py
```

三个入口共用同一个核心——Streamlit、CLI（[`interfaces/cli.py`](Model_PichiaCLM/interfaces/cli.py)）、
HTTP API（[`interfaces/api.py`](Model_PichiaCLM/interfaces/api.py)）。
ADR-0001 要求三者呈现同一份生成结果与质量判定，
所以一条候选不可能在一个入口显示合格、在另一个入口显示不合格。

## 边界

- **不预测产量。** 候选筛的是构建风险，不是表达水平。
  通过筛查不等于对实验成功的预测。
- **移植部分是上游工作** —— 见[贡献边界](#贡献边界)。
- **密码子偏好统计是描述性的**，供人比较，不作阈值。
- **候选数少于请求数是合法结果**，并且会被明确报告。

## 文档

| 文档 | 什么时候改 |
| --- | --- |
| [需求](docs/REQUIREMENTS.md) | 目标或能力边界变了 |
| [架构](docs/ARCHITECTURE.md) | 实现结构变了 |
| [执行计划](docs/EXECUTION_PLAN.md) | 进度推进了——状态的唯一权威 |
| [handoff](docs/HANDOFF.md) | 当前切片换了 |
| [ADR 索引](docs/adr/README.md) | 从不改——决策只被取代，不被改写 |

---

> 更多项目见[个人网站](https://77652189.github.io)。

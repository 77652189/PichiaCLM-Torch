# PichiaCLM - PyTorch Implementation

本项目是 **PichiaCLM** (基于深度学习的毕赤酵母密码子优化模型) 的 PyTorch 移植与优化版本。

> **🌟 来源说明 (Original Source)**
> 本项目核心架构与原始数据处理逻辑源自 [NHarini-1995/PichiaCLM](https://github.com/NHarini-1995/PichiaCLM)。原始代码基于 Keras/TensorFlow 编写，本项目将其完全重写为 PyTorch，并在推理逻辑上进行了深度的数学与生物学规则优化。

---

## 🛠️ 主要更改与优化内容 (Modifications & Improvements)

从 Keras 迁移至 PyTorch 的过程中，我们不仅复刻了原有的多任务 Seq2Seq 架构（Multi-Task Seq2Seq），还针对模型在实际推理中暴露出的问题进行了专门的修复与升级：

### 1. 框架完全迁移 (Framework Migration)
* **模型构建**：使用 PyTorch 的 `nn.Embedding` 和 `nn.GRU` 替代了原有的 Keras 网络层，并手动实现了与原版严格一致的 **点积注意力机制 (Dot-Product Attention)**。
* **损失函数**：采用了带有 `ignore_index=0` 参数的 `nn.CrossEntropyLoss`，以原生且高效的方式实现了原版 Keras 中 `mask_zero=True` 的序列填充屏蔽效果。
* **训练循环**：构建了灵活的 PyTorch 自定义训练循环，并加入了梯度裁剪 (`clip_grad_norm_`)、NaN/Inf 异常损失拦截，以及基于验证集表现的动态学习率调度器 (`ReduceLROnPlateau`)。

### 2. 修复“前缀幻觉”与对齐偏移 (Fixing Prefix Hallucination)
在早期的 Seq2Seq 测试中，我们发现模型在面对不以起始密码子（如 `M`）开头的序列时，会产生“幻觉（Hallucination）”，吞掉前几个氨基酸并输出毫不相干的密码子组合。
* **重构为无填充解码**：在推理阶段废除了固定长度的 Padding (补零)，确保双向 GRU 的隐藏状态不会受到长序列 0 的噪声污染。
* **纯步进式自回归 (Step-by-Step Autoregressive Decoding)**：将原先一次性送入整个序列的做法，改为严格的标准机器翻译推理范式——每次仅送入单步 Token 并更新 Hidden State，使得 Attention 机制能够在真实序列长度上极度聚焦。

### 3. 引入“硬性生物学规则掩码” (Hard Biological Rule Masking) —— 核心突破 🚀
普通的语言模型在翻译氨基酸时存在概率越界的风险，为了达到 **100% 的生物学翻译保真度**，我们引入了强规则掩码：
* **字典映射约束**：在解码器输出下一步密码子的概率 (Logits) 时，代码会实时截获当前位置对应的氨基酸 ID，并查询合法的密码子集合。
* **`-inf` 屏蔽机制**：对于词表中不属于该氨基酸的所有非法密码子，模型强行将其 Logits 加上负无穷大 (`-inf`)。
* **效果**：这一机制迫使模型**只能在合法的密码子选项中，挑选一个被 Attention 认为最匹配宿主（毕赤酵母）偏好性的密码子**。彻底根治了跳字、漏翻、错翻等对齐问题，实现了 AI 预测与生物学规则的完美融合。

---

## 🚀 快速开始 (Quick Start)

### 依赖环境
* Python 3.8+
* PyTorch 2.0+
* Pandas, NumPy

### 核心文件说明
* `/Model_PichiaCLM/Training/AllData/`:训练数据准备
* `Train_Test_Splitting-Torch.ipynb`: 数据清洗、统一词表构建 (SOS, EOS, UNK, PAD) 及 Tensor 序列化代码。
* `Training-Torch.ipynb`: 包含 PyTorch 模型架构、多任务训练循环 (Multi-task Training) 以及带有硬规则 Masking 的最终优化版推理脚本。
* `./Model_PichiaCLM/Training/PichiaData/2Target_AllData/Arch1-0404.weights.pt`: 训练后的权重：Val Accuracy: CDS = 88.39% | AA = 99.99%

### 推理演示代码片段
模型能够根据氨基酸序列，结合毕赤酵母表达偏好，精准生成密码子：

```python
# 初始化并加载模型
model = MultiTaskSeq2Seq().to(device)
model.load_state_dict(torch.load("Arch1-0404.weights.pt"))

# 输入目标氨基酸序列
input_aa_string = "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADES"

# 执行带有生物学约束的解码
predicted_indices = translate_aa_to_cds(model, aa_tensor, aa_id_to_codon_ids, device)
```
### 额外说明
根据训练后得到的权重推理得到的CDS序列与原论文中的数据不一致。

## 致谢 (Acknowledgments)
感谢原作者 NHarini-1995 提供的出色模型设计与毕赤酵母数据集。本 PyTorch 版本的旨在为社区提供更易于调试、扩展及二次开发的现代化深度学习实现。
Coded with Gemini3.1Pro -20260406
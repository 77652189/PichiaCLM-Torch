# PichiaCLM-Torch

中文说明 | [English](README.md)

PichiaCLM-Torch 是一个用于 *Pichia pastoris*（毕赤酵母）密码子优化的 PyTorch 实现和可部署工具。

核心流程是：

```text
蛋白质氨基酸序列 -> 适合毕赤酵母表达偏好的 CDS/DNA 序列
```

本项目不是设计新蛋白，而是把输入的氨基酸序列转换为同义编码 DNA 序列，并提供后续克隆、合成和质检所需的序列质量分析。

## 功能概览

- 单条氨基酸序列预测。
- 多条蛋白或突变体的 FASTA 批量预测。
- Streamlit 网页、FastAPI 服务和 CLI 命令行。
- 优化后 CDS 的翻译一致性检查。
- GC 含量和局部 GC 含量分析。
- CAI 和密码子使用统计，可同时对比项目训练数据和公开 Kazusa *Pichia pastoris* 密码子表。
- 连续稀有密码子、同聚碱基、重复序列、不想要的 motif、限制性酶切位点检查。
- 对外部软件二次优化后的 CDS 做质量警告，不重新预测。
- 信号肽 + 成熟蛋白拼接辅助：
  - 整体优化
  - 分段优化
- 对部分风险位点进行保守同义密码子后处理。

## 项目结构

```text
Model_PichiaCLM/
  core/
    predictor.py      # 模型加载和 AA-to-CDS 预测
    analysis.py       # 序列质量分析
    biology.py        # 密码子表、DNA 标准化、翻译检查
    fasta.py          # FASTA 解析和输出
    restriction.py    # 限制性酶切位点扫描
    fusion.py         # 信号肽 + 成熟蛋白拼接对比
    postprocess.py    # 保守同义替换后处理
  interfaces/
    cli.py            # 命令行入口
    api.py            # FastAPI 入口
    streamlit_app.py  # Streamlit 网页入口
```

设计原则：

```text
core 负责模型推理和生物序列分析
interfaces 只负责 CLI / API / Streamlit 外壳
```

## 安装

只安装核心推理依赖：

```powershell
pip install -r requirements-core.txt
```

安装 FastAPI 服务依赖：

```powershell
pip install -r requirements-api.txt
```

安装 Streamlit 页面依赖：

```powershell
pip install -r requirements-streamlit.txt
```

安装完整部署依赖：

```powershell
pip install -r requirements-deploy.txt
```

## 快速开始

### CLI

单条预测：

```powershell
python -m Model_PichiaCLM.interfaces.cli --aa MSTNPKPQR --json
```

当前内置权重的示例输出 CDS：

```text
ATGTCCACAAATCCCAAACCACAGAGA
```

FASTA 批量预测：

```powershell
python -m Model_PichiaCLM.interfaces.cli `
  --aa-fasta input.fasta `
  --analysis `
  --out-fasta output_cds.fasta `
  --out-csv report.csv
```

分析外部二次优化后的 CDS，不重新运行模型预测：

```powershell
python -m Model_PichiaCLM.interfaces.cli `
  --cds ATGTCCACAAATCCCAAACCACAGAGA `
  --expected-aa MSTNPKPQR `
  --analysis
```

### FastAPI

启动服务：

```powershell
uvicorn Model_PichiaCLM.interfaces.api:app --host 0.0.0.0 --port 8000
```

单条预测：

```powershell
Invoke-RestMethod -Method Post `
  -Uri http://127.0.0.1:8000/predict `
  -ContentType application/json `
  -Body '{"amino_acids":"MSTNPKPQR"}'
```

外部 CDS 质检：

```powershell
Invoke-RestMethod -Method Post `
  -Uri http://127.0.0.1:8000/analyze_cds `
  -ContentType application/json `
  -Body '{"cds":"ATGTCCACAAATCCCAAACCACAGAGA","expected_amino_acids":"MSTNPKPQR"}'
```

### Streamlit

启动页面：

```powershell
python -m streamlit run Model_PichiaCLM/interfaces/streamlit_app.py `
  --server.address 0.0.0.0 `
  --server.port 8501
```

本机访问：

```text
http://127.0.0.1:8501
```

Streamlit 页面包含：

- 单条预测
- FASTA 批量预测
- 二次优化 CDS 质检
- 信号肽拼接对比
- 保守后处理选项

## 模型权重

默认权重路径：

```text
Model_PichiaCLM/Training/PichiaData/2Target_AllData/Arch1-0404.weights.pt
```

CLI、API 和 Streamlit 都默认使用该路径，也可以在对应入口中手动覆盖。

## 质量检查内容

质量报告包含：

- CDS 长度和阅读框
- 与期望氨基酸序列的翻译一致性
- 内部终止密码子
- 非法 DNA 碱基
- 全局 GC 和 30 bp 滑窗局部 GC
- 基于训练数据和公开参考表的 CAI
- 密码子使用统计
- 连续稀有密码子
- 同聚碱基和重复序列
- 默认和自定义限制性酶切位点
- 用户自定义不想要的 motif

默认 GC 阈值：

```text
全局 GC: 35%-65%
局部 GC: 30 bp 窗口，25%-75%
```

默认限制性酶：

```text
EcoRI, XhoI, NotI, BamHI, HindIII, NdeI, NcoI, KpnI, XbaI, SpeI
```

## 使用建议

PichiaCLM 输出应视为候选优化 CDS。正式合成或实验使用前，建议继续复核：

- 翻译是否完全一致
- 是否符合载体和克隆策略
- 是否符合合成公司规则
- 是否满足项目特定的生物学要求

## 致谢

本仓库基于 PichiaCLM-Torch 的 PyTorch 移植和部署改造工作。PichiaCLM 的原始思想和数据处理脉络来自既有 PichiaCLM 研究代码与数据集。

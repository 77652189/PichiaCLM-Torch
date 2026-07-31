# PichiaCLM Handoff

## 当前切片

完成文档治理基线及 ADR-0001：最终候选的保守相对基准筛选是 `core` 的唯一规则；CAI 只作为双参考审查指标，不是单独的合格阈值。

## 下一步

等待用户选择并授权一个具体产品、科学或部署切片。收到授权后，先将其进度与停止条件写入 `docs/EXECUTION_PLAN.md`，再执行最小可验证改动。

## 必读材料

1. `docs/REQUIREMENTS.md`
2. `docs/ARCHITECTURE.md`
3. `docs/adr/ADR-0001-qualified-candidate-acceptance.md`
4. `docs/EXECUTION_PLAN.md`

## 验证方式

运行 `python -m pytest -q tests/test_core_features.py tests/test_docs_governance.py`。若涉及接口、模型或部署，再增加与该切片对应的端到端验证。

## 硬约束

- PichiaCLM 输出是候选 CDS 和审查信息，不是表达产量、分泌效率或湿实验成功的预测。
- 不修改训练数据或模型权重，除非用户明确授权一个包含数据来源、评估门控和回滚方式的切片。
- 最终候选不得有关键问题、不得增加风险警告，也不得增加可避免的最低偏好密码子；CAI 不单独决定合格状态。
- 不改变远端可见性、不提交、不推送，除非用户明确要求。

```yaml
current_slice: documentation_governance
slice_status: done
authorization_status: awaiting_user
verification_status: passed
```

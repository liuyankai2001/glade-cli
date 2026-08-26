# RetroPath P11.1 隐藏 KEGG 反应基线报告

- Benchmark：`p11_1_hidden_kegg_pilot12`
- Run ID：`20260825T090033Z_4f69c442`
- 创建时间：`2026-08-25T09:00:33.556515+00:00`
- 数据集 SHA-256：`2d21764cdb4652231fef5a655ab8207ae35ad981b688910eee0283b7ddac51d9`
- Git commit：`3ea3b3f56944a5e19e3139480b4364497f791e0e`
- 工作树状态：`dirty`

> 本报告是 P11.1 正向恢复能力测试。金标准反应的来源规则仍保留在 RR02 中，
> 因此结果表示可恢复上限，不是未知反应假阳性率；来源规则排除属于 P11.2。

## 执行概况

| Profile | 计划任务 | 可评测任务 | 完成率 | Scope 命中率 | 截断率 | 运行时间中位数 |
|---|---:|---:|---:|---:|---:|---:|
| controlled | 12 | 7 | 58.3% | 42.9% | 0.0% | 44.6s |
| full_a0 | 12 | 6 | 50.0% | 50.0% | 0.0% | 41.2s |

## 核心恢复率

下表使用可评测案例作为分母；CSV 同时保留 `all_selected` 保守下限。

| Profile | 指标 | @1 | @3 | @5 | @10 |
|---|---|---:|---:|---:|---:|
| controlled | 来源模板恢复 | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) |
| controlled | 平衡计量恢复 | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) |
| controlled | 严格 GEM 恢复 | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) |
| controlled | 精确正式 solution | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) | 14.3% (1/7) |
| full_a0 | 来源模板恢复 | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) |
| full_a0 | 平衡计量恢复 | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) |
| full_a0 | 严格 GEM 恢复 | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) |
| full_a0 | 精确正式 solution | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) | 16.7% (1/6) |

## 首次命中排名质量

| Profile | 模板 MRR | 平衡 MRR | GEM MRR | 精确 solution MRR |
|---|---:|---:|---:|---:|
| controlled | 0.143 | 0.143 | 0.143 | 0.143 |
| full_a0 | 0.167 | 0.167 | 0.167 | 0.167 |

## 双 Profile 配对比较

- 可比较案例：6
- 两种 profile 均精确恢复：1
- 仅 controlled 恢复：0
- 仅 full_a0 恢复：0
- 两种 profile 均未恢复：5

## 单案例结果

| Case | EC 类 | Profile | 状态 | Scope | 模板排名 | 平衡排名 | GEM 排名 | 精确排名 | 耗时 |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|
| ec1_r01299 | 1 | controlled | completed | 是 | - | - | - | - | 6.8s |
| ec1_r01299 | 1 | full_a0 | completed | 是 | - | - | - | - | 10.5s |
| ec1_r01574 | 1 | controlled | failed | 否 | - | - | - | - | 3.7s |
| ec1_r01574 | 1 | full_a0 | failed | 否 | - | - | - | - | 4.6s |
| ec2_r00006 | 2 | controlled | completed | 是 | - | - | - | - | 76.2s |
| ec2_r00006 | 2 | full_a0 | failed | 否 | - | - | - | - | 73.4s |
| ec2_r10062 | 2 | controlled | failed | 否 | - | - | - | - | 0.2s |
| ec2_r10062 | 2 | full_a0 | failed | 否 | - | - | - | - | 1.6s |
| ec3_r00913 | 3 | controlled | completed | 是 | 1 | 1 | 1 | 1 | 49.7s |
| ec3_r00913 | 3 | full_a0 | completed | 是 | 1 | 1 | 1 | 1 | 10.9s |
| ec3_r10629 | 3 | controlled | completed | 否 | - | - | - | - | 44.2s |
| ec3_r10629 | 3 | full_a0 | completed | 否 | - | - | - | - | 41.8s |
| ec4_r02305 | 4 | controlled | completed | 否 | - | - | - | - | 39.6s |
| ec4_r02305 | 4 | full_a0 | completed | 否 | - | - | - | - | 42.3s |
| ec4_r10597 | 4 | controlled | completed | 否 | - | - | - | - | 80.8s |
| ec4_r10597 | 4 | full_a0 | completed | 是 | - | - | - | - | 65.8s |
| ec5_r01057 | 5 | controlled | failed | 否 | - | - | - | - | 105.8s |
| ec5_r01057 | 5 | full_a0 | failed | 否 | - | - | - | - | 74.9s |
| ec5_r02007 | 5 | controlled | completed | 否 | - | - | - | - | 44.6s |
| ec5_r02007 | 5 | full_a0 | completed | 否 | - | - | - | - | 40.6s |
| ec6_r01322 | 6 | controlled | failed | 否 | - | - | - | - | 0.1s |
| ec6_r01322 | 6 | full_a0 | skipped | 否 | - | - | - | - | 0.0s |
| ec6_r01991 | 6 | controlled | failed | 否 | - | - | - | - | 44.1s |
| ec6_r01991 | 6 | full_a0 | failed | 否 | - | - | - | - | 42.9s |

## 失败与限制

- `ec1_r01574/controlled`：pipeline_error / retropath_parse_failed — artifact_inconsistent: scope has no iteration-0 transformation rooted at the P2 target
- `ec1_r01574/full_a0`：pipeline_error / retropath_parse_failed — artifact_inconsistent: scope has no iteration-0 transformation rooted at the P2 target
- `ec2_r00006/full_a0`：pipeline_error / retropath_parse_failed — sink_identity_mismatch: service reported source_in_sink but target is not in the P2 sink
- `ec2_r10062/controlled`：infrastructure_error / retropath_timeout — RetroPath execution exceeded 3600 seconds
- `ec2_r10062/full_a0`：infrastructure_error / retropath_execution_failed — service restarted while the job was running
- `ec5_r01057/controlled`：pipeline_error / retropath_parse_failed — artifact_inconsistent: scope has no iteration-0 transformation rooted at the P2 target
- `ec5_r01057/full_a0`：pipeline_error / retropath_parse_failed — sink_identity_mismatch: service reported source_in_sink but target is not in the P2 sink
- `ec6_r01322/controlled`：infrastructure_error / retropath_service_unavailable — submission_uncertain: a previous submission may have reached the service; use force=True only after checking the service
- `ec6_r01322/full_a0`：infrastructure_circuit_breaker / controlled_profile_infrastructure_failure — full_a0 was not submitted because the paired controlled profile ended with infrastructure error retropath_service_unavailable
- `ec6_r01991/controlled`：pipeline_error / retropath_parse_failed — artifact_inconsistent: scope has no iteration-0 transformation rooted at the P2 target
- `ec6_r01991/full_a0`：pipeline_error / retropath_parse_failed — artifact_inconsistent: scope has no iteration-0 transformation rooted at the P2 target

本试运行只有 12 个反应案例。比例使用 Wilson 95% 区间，运行时间使用中位数和 IQR；
不进行显著性检验，也不基于本报告自动修改生产参数。

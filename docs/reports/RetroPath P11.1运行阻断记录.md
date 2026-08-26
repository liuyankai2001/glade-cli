# RetroPath P11.1 真实回测运行与恢复记录

## 最终状态

- Run ID：`20260825T090033Z_4f69c442`
- 运行目录：`tests/.retropath_benchmark_runtime/20260825T090033Z_4f69c442`
- 计划任务：12 个案例 × 2 个 profile，共 24 个核心任务
- 最终记录：24 个终态任务；13 个可评测、10 个失败、1 个熔断跳过
- 正式报告：`docs/reports/RetroPath P11.1基线报告.md`
- Run 状态：`completed_with_failures`

## 故障与恢复经过

首次真实回测中，单 worker 的 KNIME 长作业超过一小时未释放，Docker Desktop Linux
Engine 管理 API 持续返回 HTTP 500。用户重启 Docker Desktop 后，原 run 从断点继续，
前 20 个任务得到终态记录。第 21 个 `ec6_r01322/controlled` 再次令 Linux Engine
端点退出；Docker Desktop 日志持续报告无法连接 `192.168.65.7:2376`。

为避免重启后重复提交同一危险作业，评测器增加了两层保护：

1. 哈希有效的 `completed`、`failed`、`skipped` 默认都作为终态恢复；仅显式
   `--retry-failed` 才重试失败任务；
2. controlled profile 发生基础设施故障时，同案例 full-A0 记录为熔断跳过，不再
   立即提交相同目标。

随后使用 Docker Desktop 官方 CLI 恢复 Engine。第 21 个作业按“已提交、终态不确定”
记录为基础设施失败，第 22 个任务熔断；第 23、24 个任务正常提交并返回，最终报告器
校验全部 task/artifact SHA-256 后生成正式报告。

## 基线摘要

| 指标 | controlled | full-A0 |
|---|---:|---:|
| 计划任务 | 12 | 12 |
| 可评测任务 | 7 | 6 |
| 运行完成率 | 58.3% | 50.0% |
| 精确恢复（all-selected） | 1/12（8.3%） | 1/12（8.3%） |
| 精确恢复（evaluable） | 1/7（14.3%） | 1/6（16.7%） |
| 可评测任务运行时间中位数 | 44.6 s | 41.2 s |

失败漏斗包含 7 个 `retropath_parse_failed`、3 个基础设施失败和 1 个配对熔断。
P4 失败主要是缺少以 P2 target 为根的 iteration-0 transformation，或服务把结果报告为
`source_in_sink`、但 target 并不在提交的 sink 中。正式报告同时保留 all-selected
保守分母和 evaluable 分母；失败没有从保守恢复率中消失。

## 后续建议

P11.2 之前先处理两类问题：

1. 明确 RetroPath `results.csv` 中 source/sink 与 iteration 的实际语义，修复或兼容
   P4 artifact 一致性判断；
2. 为 Docker/KNIME worker 增加真正的进程级资源上限与硬超时，避免单一化合物令整个
   Linux Engine 退出。

不要直接用本次 14.3%/16.7% 的 evaluable 恢复率调低假阳性阈值；当前主要瓶颈仍是
可运行性和原始网络解释，而不是候选排序阈值。

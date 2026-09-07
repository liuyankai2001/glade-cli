# RetroPath 本地服务使用说明

## 1. 运行栈

本地服务把 RetroPath 的非 Python 运行环境与 GLADE 的 uv 环境隔离：

| 组件 | 固定版本 |
|---|---|
| retropath2_wrapper | 3.9.1 |
| RetroPath2.0 KNIME workflow | r20260212 |
| KNIME | 4.7.0 |
| KNIME RDKit Nodes | 4.9.1 |
| KNIME Chromium compatibility | isolated OpenSSL 1.0.2k prefix |
| RetroRules | rr02-rp2-hs |
| RR02 retro SHA-256 | e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f |

服务只监听 `127.0.0.1:8765`，不会暴露到局域网。GLADE 暂未在 P0 阶段调用该接口。

上游 3.9.1 conda 包中的 Python 元数据会误报为 3.9.0；健康检查以
conda 安装清单确认实际版本，并额外返回 `wrapper_reported_version` 保留该差异。

## 2. 前置检查

确认 Docker Desktop 已启动：

```powershell
docker version
```

确认 RR02 逆向规则存在：

```powershell
Get-Item data\retropath\rules\rr02\retrorules_rr02_rp2_flat_retro.csv
```

## 3. 构建与启动

首次构建会下载 KNIME、化学节点和 RDKit 节点，耗时和镜像体积都明显高于普通 Python 服务：

```powershell
docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d retropath
```

查看状态和日志：

```powershell
docker compose -f compose.retropath.yml ps
docker compose -f compose.retropath.yml logs -f retropath
```

健康检查：

```powershell
Invoke-RestMethod http://127.0.0.1:8765/health
```

只有 `ready` 为 `true` 时才能提交任务。

## 4. 提交任务

source 和 sink 均为 UTF-8 CSV，前两列必须是 `Name,InChI`。source 必须恰好包含一个化合物。

```powershell
curl.exe -X POST http://127.0.0.1:8765/v1/jobs `
  -F "source_file=@source.csv;type=text/csv" `
  -F "sink_file=@sink.csv;type=text/csv" `
  -F "max_steps=3" `
  -F "topx=100" `
  -F "dmin=2" `
  -F "dmax=16"
```

返回结果中的 `job_id` 用于查询：

```powershell
Invoke-RestMethod http://127.0.0.1:8765/v1/jobs/<job_id>
Invoke-RestMethod http://127.0.0.1:8765/v1/jobs/<job_id>/results
```

状态包括：`queued`、`running`、`succeeded`、`no_solution`、`source_in_sink`、`failed` 和 `timed_out`。

## 5. 测试

运行容器内单元测试：

```powershell
docker compose -f compose.retropath.yml run --rm retropath pytest -q /opt/service/tests
```

服务启动后运行快速接口测试：

```powershell
$env:RETROPATH_SERVICE_URL = "http://127.0.0.1:8765"
uv run pytest -q tests/test_retropath_local_service.py
```

真实 KNIME/RR02 冒烟测试：

```powershell
$env:RETROPATH_FUNCTIONAL = "1"
uv run pytest -q tests/test_retropath_local_service.py
```

## 6. 停止

```powershell
docker compose -f compose.retropath.yml down
```

不要默认使用 `down -v`，否则会删除保存任务和结果的命名卷。

## 7. 内存超限排查

出现 `RetroPath cgroup working set exceeded ... for 3 consecutive samples`
时，服务监控器检测到工作集连续超出软上限，主动终止任务。
服务端 `failure_code` 为 `resource_exhausted`，流水线状态仍为
`failed`。这不能解释为“没有合成路线”。GLADE 会检查服务已列入清单并完成
哈希校验的 `results.csv`，只恢复其中从目标到可信 sink 的完整闭合路线。

若恢复到候选，`pipeline_result.json` 使用
`retropath_partial_candidates_found`，并记录 `search_complete=false` 和
`interrupted_result_recovered=true`。这些候选可以继续接受严格 GEM 验证，
但不能说明搜索空间已经穷尽。若没有闭合路线、结果损坏或失败原因不是内存超限/
墙钟超时，流水线仍返回原来的执行失败。

查看本次输出目录的 `raw/service_run_manifest.json`：

- `parameters`：实际提交的 `max_steps`、`topx` 等参数。
- `resource_telemetry`：工作集峰值、软上限、连续超限次数和运行时间。
- `failure_code`、`return_code`：失败原因和进程退出码。

`pipeline_result.json` 同时保留 `service_failure`（任务 ID、失败码、退出码、
提交参数、服务监控记录路径）。没有恢复出候选时，内存超限错误还提供
`recovery_hint`。

例如 C05432 的一次运行使用 `max_steps=5`、`topx=100`，约 112 秒后工作集
峰值达到 10.35 GiB，触发 10 GiB 软上限。监控已减去可回收的
`inactive_file`，该次原始内存与工作集仅相差约 5.5 MiB。

如果可以接受缩小搜索范围，可在原命令中将 `--step 5` 改为 `--step 3`，
其余参数保持一致。`--step` 控制 RetroPath 逆合成步数，`-d` 控制底盘扩展
深度；两者不同。降低步数不保证解决内存问题，也不能证明更长路线不存在。
`--input` 仍只传 `inputs` 目录下的文件名。

如果必须保留搜索范围，先检查实际资源配置：

```powershell
docker info --format '{{.MemTotal}} {{.NCPU}}'
docker compose -f compose.retropath.yml exec -T retropath cat /sys/fs/cgroup/memory.max
docker compose -f compose.retropath.yml exec -T retropath cat /sys/fs/cgroup/cpu.max
docker compose -f compose.retropath.yml exec -T retropath cat /opt/knime/knime_4.7.0/knime.ini
```

当前镜像设置 `-Xmx2048m`，它只限制 Java 堆；cgroup 工作集还包含原生库、
线程和其他进程等内存。若日志显示大量 KNIME 工作线程，应进一步排查运行
并发对峰值内存的影响。不要仅根据 Java 堆上限判断容器内存需求。

软上限由 `.env` 中的 `RETROPATH_MEMORY_LIMIT_BYTES` 控制，容器硬上限由
`RETROPATH_CONTAINER_MEMORY_LIMIT` 控制；调整和重建容器的方法见
[GLADE 部署教程](GLADE部署教程.md)。必须同时核对 Docker/WSL 总内存和
宿主机余量。例如硬上限已有 11 GiB，而 Docker 总内存约 11.68 GiB 时，
继续提高软上限的空间很小，应优先评估搜索规模及并发。

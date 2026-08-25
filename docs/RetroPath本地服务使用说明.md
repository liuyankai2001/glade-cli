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

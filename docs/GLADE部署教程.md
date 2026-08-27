# GLADE 详细部署教程

本文档说明如何在一台新的 Windows 计算机上从零部署 GLADE。教程以 Windows 10/11、
PowerShell 7、Python 3.12、`uv` 和 Docker Desktop 为主；所有命令都应在 GLADE 项目
根目录执行。

部署完成后的实际使用方法见 [GLADE 用户使用说明](GLADE用户使用说明.md)。

## 1. 先选择需要部署的功能

GLADE 的不同功能依赖不同组件，不需要为了运行基础分析而一次安装所有服务。

| 部署级别 | 可以使用的功能 | 需要的组件 |
|---|---|---|
| 基础路线分析 | `chassis`、`expand`、KEGG `gap`、KEGG 路线验证 | Python 环境、GEM、培养基、KEGG 等数据库网络 |
| RetroPath 搜索 | 基础功能加 `gap --retropath` | Docker Desktop、RetroPath 规则文件、本地 RetroPath 服务 |
| RetroPath 路线验证 | RetroPath 搜索加预测路线 GEM 验证 | MNXref v3.0 本地索引 |
| 完整设计流程 | 主酶、CDS、表达元件、质粒和组装 | CodonTransformer 模型、远端 Milvus；需要结构回退检索时使用 Selenzyme，研究功能还需要大模型接口 |

建议先完成基础部署和冒烟测试，再按需要安装 RetroPath、MNXref 和完整设计流程所需的
组件。这样出现问题时更容易判断是 Python 环境、数据文件还是外部服务造成的。

## 2. 部署完成标准

基础部署至少应满足：

- `uv run python main.py -h` 能显示 GLADE 命令；
- Python 依赖检查通过；
- `data/gem_models/iML1515.json` 和 `data/mediums/default_medium.json` 存在；
- 能够对测试配置运行 `chassis` 并查看结果。

如果部署 RetroPath，还应满足：

- Docker Compose 配置检查通过；
- `http://127.0.0.1:8765/health` 返回 `ready: true`；
- 健康检查中的规则文件 SHA-256 与本教程一致。

如果部署完整流程，还应满足：

- CodonTransformer 本地模型文件完整；
- MNXref 状态检查通过；
- Milvus 中存在 `expression_parts_v3` 和 `plasmid_templates_v2`；
- 所需的大模型和 Selenzyme 地址可以从部署机器访问。

## 3. 系统和资源要求

### 3.1 操作系统与基础软件

推荐环境：

- 64 位 Windows 10 或 Windows 11；
- PowerShell 7；
- Git；
- Python 3.12，项目不支持 Python 3.13；
- `uv`；
- 需要 RetroPath 时安装 Docker Desktop，并使用 Linux 容器。

官方安装文档：

- [uv 安装说明](https://docs.astral.sh/uv/getting-started/installation/)
- [Docker Desktop for Windows](https://docs.docker.com/desktop/setup/install/windows-install/)
- [Git for Windows](https://git-scm.com/download/win)

### 3.2 内存和磁盘

基础路线分析和 CDS 生成可以只使用 CPU，CUDA 不是必需条件。

RetroPath 容器在 `compose.retropath.yml` 中设置了 7 GB 内存上限，并会在内存持续超过
6 GiB 时中止任务。因此 Docker Desktop 必须能获得至少 7 GB 内存；宿主机还需要为
Windows、Python 和文件缓存保留余量。实际部署时建议使用 12–16 GiB 或更多物理内存。

以下内容不会随普通 Git 克隆完整获得，并且会占用额外磁盘空间：

- RetroPath 规则文件：当前所需逆向规则文件约 204 MB，完整压缩包解压后更大；
- MNXref 本地索引：约 80 MB，安装时还需要下载和处理更大的源文件；
- CodonTransformer 模型：约 342 MB；
- RetroPath Docker 镜像：包含 KNIME、化学节点和 Conda 环境，需要数 GB 空间。

`model/bge-m3` 当前不参与 GLADE 的表达元件或质粒推荐，部署当前流程时不需要下载。

### 3.3 网络访问

按使用的功能放行相应地址：

| 用途 | 地址或服务 |
|---|---|
| 安装 Python 依赖 | Python 包索引和 `uv` 使用的软件源 |
| KEGG 路线与结构 | `https://rest.kegg.jp` |
| 反应和酶证据 | Rhea、UniProt |
| RetroPath Docker 构建 | Docker Registry、Conda Forge、KNIME 更新站点、Zenodo |
| MNXref 安装 | `https://www.metanetx.org` |
| CodonTransformer 下载 | `https://huggingface.co` |
| 表达元件和质粒 | 部署方提供的 Milvus 地址 |
| 文献和辅助蛋白研究 | `.env` 中配置的大模型接口 |
| 结构相似酶检索 | `.env` 中配置的 Selenzyme 地址 |

本地 RetroPath 服务只绑定 `127.0.0.1:8765`，不需要也不应暴露到局域网或公网。

## 4. 获取项目源码

如果尚未获得项目目录，使用部署方提供的仓库地址：

```powershell
$gladeRepositoryUrl = '<GLADE仓库地址>'
git clone $gladeRepositoryUrl glade
Set-Location .\glade
```

如果已经拿到完整项目目录，直接进入根目录。以下文件应位于当前目录：

```powershell
Get-Item main.py, pyproject.toml, uv.lock, compose.retropath.yml
```

后续命令不要在 `src`、`docs` 或其他子目录中执行，否则相对路径和 `.env` 位置会不一致。

## 5. 安装 uv 和 Python 3.12

### 5.1 安装 uv

若尚未安装 `uv`，可以使用官方 PowerShell 安装脚本：

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

关闭并重新打开 PowerShell，然后确认：

```powershell
uv --version
```

如果单位安全策略不允许执行在线脚本，请按照官方文档下载并安装 `uv`，不要绕过单位的
安全策略。

### 5.2 安装项目要求的 Python

```powershell
uv python install 3.12
uv python find 3.12
```

项目在 `pyproject.toml` 中要求 `>=3.12,<3.13`。不要使用 3.11 或 3.13 创建环境。

### 5.3 同步锁定依赖

```powershell
uv sync --frozen
```

`--frozen` 会直接使用仓库中的 `uv.lock`，避免部署时意外更新依赖版本。成功后项目根目录
会出现 `.venv`。

检查 Python 和依赖：

```powershell
uv run python --version
uv pip check
uv run python main.py -h

@'
import Bio
import cobra
import CodonTransformer
import dnachisel
import httpx
import numpy
import pandas
import pymilvus
import rdkit
import RNA
print("GLADE Python dependencies: OK")
'@ | uv run python -
```

后续推荐统一使用 `uv run python ...`，不需要手动激活虚拟环境。如果确实需要激活：

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
& .\.venv\Scripts\Activate.ps1
```

激活后可以把 `uv run python` 简写为 `python`。

## 6. 检查基础数据和目录

基础分析依赖默认大肠杆菌 GEM 和培养基：

```powershell
$baseFiles = @(
  'data\gem_models\iML1515.json',
  'data\mediums\default_medium.json'
)
$missingBaseFiles = $baseFiles | Where-Object { -not (Test-Path -LiteralPath $_) }
if ($missingBaseFiles) {
  throw "缺少基础数据文件: $($missingBaseFiles -join ', ')"
}
```

确认运行目录存在：

```powershell
New-Item -ItemType Directory -Force -Path inputs, outputs, cache, model | Out-Null
```

这些目录的用途是：

| 目录 | 用途 |
|---|---|
| `inputs` | 用户输入配置和手动上传的蛋白/CDS 文件 |
| `outputs` | 每个目标化合物的全部项目结果 |
| `cache` | KEGG、Rhea、UniProt、结构和文献等缓存 |
| `model` | CodonTransformer 等本地模型 |
| `data` | GEM、培养基、RetroPath 规则和 MNXref 数据 |

`inputs`、`outputs`、`cache`、`model` 和 `data/retropath` 都不会被普通 Git 提交自动备份。
生产部署应单独制定备份策略。

## 7. 配置 `.env`

不是所有命令都需要 `.env`。基础的 `chassis`、`expand`、KEGG 路线搜索和 GEM 验证
不依赖大模型或 Milvus。

需要完整设计流程时，在项目根目录创建 `.env`：

```dotenv
# 辅助蛋白研究和 --literature-search 使用
MODEL_PROVIDER=openai
AGENT_LLM_MODEL=<模型名称>
API_KEY=<API密钥>
BASE_URL=<兼容OpenAI API的服务地址>

# 可选：蛋白研究缓存目录；相对路径按项目根目录解析
CACHE_DIR=cache/protein_supply

# 可选：为数据库请求提供联系信息
GLADE_CONTACT_EMAIL=<联系邮箱>

# RetroPath 预测步骤和数据库未覆盖步骤的结构相似酶检索
SELENZYME_REST_URL=<Selenzyme REST服务地址>

# 表达元件和质粒骨架推荐
MILVUS_HOST=<Milvus主机或完整HTTP地址>
MILVUS_PORT=19530
MILVUS_TOKEN=<可选令牌>
MILVUS_DB_NAME=<可选数据库名>
```

配置规则：

- `MODEL_PROVIDER` 当前只支持 `openai`；
- `BASE_URL` 必须兼容 OpenAI API；
- 如果 `MILVUS_HOST` 已经写成 `http://...` 或 `https://...`，程序直接使用该完整地址；
- 如果 `MILVUS_HOST` 只写主机名或 IP，程序会使用 `MILVUS_PORT` 拼接地址；
- 不要把真实 `.env` 提交到 Git，也不要在日志、截图或问题报告中公开 `API_KEY`、
  `MILVUS_TOKEN` 等秘密；
- 当前 GLADE 源码不读取 MySQL 配置，不需要为了运行 GLADE 添加 MySQL 变量。

### 7.1 Milvus 前置条件

当前仓库只包含 Milvus 客户端，不包含创建和填充生产集合的脚本。部署方必须提供已经
填充并通过审计的集合：

| 集合 | 用途 |
|---|---|
| `expression_parts_v3` | promoter、RBS 和 terminator 推荐 |
| `plasmid_templates_v2` | 质粒骨架推荐 |

可以用下面的只读脚本检查连接和集合是否存在：

```powershell
@'
import os
from dotenv import load_dotenv
from pymilvus import MilvusClient

load_dotenv(".env")
host = (os.getenv("MILVUS_HOST") or "").strip()
port = (os.getenv("MILVUS_PORT") or "19530").strip()
if not host:
    raise SystemExit("MILVUS_HOST 未配置")
uri = host.rstrip("/") if host.startswith(("http://", "https://")) else f"http://{host}:{port}"
kwargs = {"uri": uri}
token = (os.getenv("MILVUS_TOKEN") or "").strip()
db_name = (os.getenv("MILVUS_DB_NAME") or "").strip()
if token:
    kwargs["token"] = token
if db_name:
    kwargs["db_name"] = db_name
client = MilvusClient(**kwargs)
for name in ("expression_parts_v3", "plasmid_templates_v2"):
    print(f"{name}: {client.has_collection(name)}")
'@ | uv run python -
```

两个集合都应输出 `True`。程序运行时还会继续检查字段和数据版本，仅“集合存在”不代表
集合内容一定满足要求。

## 8. 安装 CodonTransformer 本地模型

只有执行 `protein-to-cds` 并需要把氨基酸序列优化为 CDS 时才需要该模型。用户直接上传
的 CDS 会跳过密码子优化，但路线中的其他蛋白如果仍是氨基酸序列，仍需要本地模型。

模型必须位于固定目录：

```text
model/CodonTransformer/
├── config.json
├── model.safetensors
├── tokenizer.json
├── tokenizer_config.json
└── ...
```

使用随 Python 依赖安装的 Hugging Face CLI 下载：

```powershell
uv run hf download adibvafa/CodonTransformer --local-dir model/CodonTransformer
```

模型来源：[adibvafa/CodonTransformer](https://huggingface.co/adibvafa/CodonTransformer)。
若部署方提供了已经审核的模型快照，应优先复制该快照，并记录其 revision 或文件哈希。

检查关键文件：

```powershell
@'
from pathlib import Path

model_dir = Path("model/CodonTransformer")
required = ["model.safetensors", "config.json", "tokenizer.json"]
missing = [name for name in required if not (model_dir / name).is_file()]
assert not missing, f"missing: {missing}"
print("CodonTransformer model: OK")
'@ | uv run python -
```

CPU 可以运行密码子优化，只是速度较慢。只有已经正确安装兼容 CUDA 的 PyTorch 和显卡
驱动时才使用 `--device cuda`；否则使用 `--device cpu` 或默认的 `--device auto`。

## 9. 安装 RetroPath 规则文件

这一节仅在需要 `gap --retropath` 时执行。

RetroPath 使用固定的逆向规则文件：

```text
data/retropath/rules/rr02/retrorules_rr02_rp2_flat_retro.csv
```

该目录被 Git 忽略，新克隆的仓库通常不会包含规则文件。规则来源为
[RetroRules release rr02-rp2-hs](https://doi.org/10.5281/zenodo.5828017)，许可证为
CC BY 4.0。部署到生产或商业环境前，还应确认规则引用的各来源数据库许可符合使用场景。

下载并解压官方归档：

```powershell
$retroRulesDir = Join-Path (Get-Location) 'data\retropath\rules\rr02'
$retroRulesArchive = Join-Path (Get-Location) 'cache\retrorules_rr02_rp2_hs.tar.gz'
$retroRulesExtractDir = Join-Path (Get-Location) 'cache\retrorules_rr02_extracted'
New-Item -ItemType Directory -Force -Path $retroRulesDir | Out-Null
New-Item -ItemType Directory -Force -Path (Split-Path $retroRulesArchive) | Out-Null
New-Item -ItemType Directory -Force -Path $retroRulesExtractDir | Out-Null

Invoke-WebRequest `
  -Uri 'https://zenodo.org/api/records/5828017/files/retrorules_rr02_rp2_hs.tar.gz/content' `
  -OutFile $retroRulesArchive

tar -xzf $retroRulesArchive -C $retroRulesExtractDir

$retroRulesSource = Get-ChildItem -LiteralPath $retroRulesExtractDir -Recurse -File `
  -Filter 'retrorules_rr02_rp2_flat_retro.csv' | Select-Object -First 1
if ($null -eq $retroRulesSource) {
  throw '下载的归档中没有找到 RetroPath 逆向规则文件'
}
Copy-Item -LiteralPath $retroRulesSource.FullName `
  -Destination (Join-Path $retroRulesDir 'retrorules_rr02_rp2_flat_retro.csv') `
  -Force
```

检查文件和 SHA-256：

```powershell
$retroRulesFile = Join-Path $retroRulesDir 'retrorules_rr02_rp2_flat_retro.csv'
$expectedRetroRulesHash = 'e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f'
if (-not (Test-Path -LiteralPath $retroRulesFile)) {
  throw "没有找到 RetroPath 逆向规则文件: $retroRulesFile"
}
$actualRetroRulesHash = (Get-FileHash -LiteralPath $retroRulesFile -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actualRetroRulesHash -ne $expectedRetroRulesHash) {
  throw "RetroPath 规则文件 SHA-256 不匹配: $actualRetroRulesHash"
}
Write-Host 'RetroPath rules: OK'
```

哈希不匹配时不要修改 `compose.retropath.yml` 绕过检查，应重新从指定来源下载正确文件。

## 10. 构建并启动 RetroPath 本地服务

### 10.1 检查 Docker

启动 Docker Desktop，确认使用 Linux 容器，然后执行：

```powershell
docker version
docker compose version
docker compose -f compose.retropath.yml config --quiet
```

第三条命令没有输出并返回退出码 0，表示 Compose 文件语法有效。

### 10.2 首次构建

```powershell
docker compose -f compose.retropath.yml build retropath
```

首次构建会下载基础镜像、Conda 软件包、KNIME 4.7、化学节点和 RDKit 节点，明显慢于
普通 Python 镜像。构建过程中需要持续访问 Docker Registry、Conda Forge、KNIME 和
Zenodo。网络中断时可以在网络恢复后重新执行同一命令，Docker 会复用已完成的层。

### 10.3 启动和健康检查

```powershell
docker compose -f compose.retropath.yml up -d retropath
docker compose -f compose.retropath.yml ps
```

容器启动后执行：

```powershell
$retroPathHealth = Invoke-RestMethod http://127.0.0.1:8765/health
$retroPathHealth | ConvertTo-Json -Depth 5
if (-not $retroPathHealth.ready) {
  throw "RetroPath 服务尚未就绪: $($retroPathHealth.errors -join '; ')"
}
```

成功结果至少应满足：

- `ready` 为 `true`；
- `wrapper_version` 为 `3.9.1`；
- `workflow_version` 为 `r20260212`；
- `knime_version` 为 `4.7.0`；
- `rules_sha256` 为
  `e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f`。

查看最近日志：

```powershell
docker compose -f compose.retropath.yml logs --tail 100 retropath
```

### 10.4 日常启动和停止

镜像已经构建后，日常只需：

```powershell
docker compose -f compose.retropath.yml up -d retropath
```

使用结束后停止：

```powershell
docker compose -f compose.retropath.yml down
```

不要默认执行 `down -v`，否则会删除保存 RetroPath 任务和结果的 Docker 命名卷。

## 11. 安装 MNXref 验证索引

这一节仅在需要验证 RetroPath 预测路线时执行。单纯搜索、查看和写入 RetroPath 路线
不要求先安装 MNXref。

规则文件安装完成后执行：

```powershell
uv run python -m src.pathway_analyze.retropath_mnxref install
```

安装器会从 MetaNetX 下载 MNXref v3.0 源文件，校验官方 MD5，再构建只包含当前规则
所需反应、化合物和映射的 SQLite 索引。下载和建库会耗时，并需要足够临时磁盘空间。

检查安装：

```powershell
uv run python -m src.pathway_analyze.retropath_mnxref status
```

成功时输出中应包含：

- `mnxref_version` 为 `3.0`；
- `rr02_sha256` 与第 9 节规则文件哈希一致；
- `index_path` 指向
  `data/retropath/mnxref/3.0/mnxref_rr02_subset.sqlite3`；
- `counts` 中存在已索引的反应、化合物和规则映射数量。

不要手工编辑 SQLite 或 manifest。规则文件变化后应重新运行安装器，使索引与规则重新
绑定。

## 12. 创建输入并运行基础冒烟测试

### 12.1 创建测试配置

下面的 PowerShell 7 命令会创建 `inputs/deploy_smoke.json`：

```powershell
@'
{
  "target_name": "C00811"
}
'@ | Set-Content -LiteralPath inputs\deploy_smoke.json -Encoding utf8NoBOM
```

所有 GLADE 命令的 `--input` 后只写 `inputs` 目录下的文件名。正确写法是
`-i deploy_smoke.json`，不要写完整路径或 `inputs/deploy_smoke.json`。

### 12.2 运行底盘分析

```powershell
uv run python main.py chassis -i deploy_smoke.json
uv run python main.py info -i deploy_smoke.json --chassis
```

成功后应生成：

```text
outputs/C00811/chassis_result/
├── producible_kegg_compounds.csv
└── analyze_chassis_metabolites_summary.csv
```

这一步不需要 Docker、MNXref、Milvus、Selenzyme 或大模型接口。如果这一步失败，应先
修复 Python 环境、GEM 或培养基问题，不要继续排查 RetroPath。

### 12.3 可选：验证 KEGG 网络访问

```powershell
uv run python main.py gap -i deploy_smoke.json -d 0
uv run python main.py info -i deploy_smoke.json --gap -d 0
```

该步骤会访问 KEGG。能正常完成或明确报告“没有候选路线”，都说明命令已运行；网络
超时、DNS 或 TLS 错误则表示外部网络尚未配置好。

### 12.4 可选：验证 RetroPath 完整调用

确认本地服务 `ready: true` 后执行：

```powershell
uv run python main.py gap -i deploy_smoke.json --retropath -d 0
uv run python main.py info -i deploy_smoke.json --retropath -d 0
```

真实 RetroPath 任务可能运行较长时间，也可能合法地返回没有完整候选。部署验收应关注
服务是否被成功调用、任务是否得到明确终态，以及 `info --retropath` 是否能解释结果，
不应把“必须找到候选路线”作为环境是否部署成功的唯一判断标准。

## 13. 完整部署验收清单

### 13.1 基础环境

- [ ] 当前目录包含 `main.py`、`pyproject.toml` 和 `uv.lock`；
- [ ] `uv --version` 可用；
- [ ] `uv run python --version` 显示 Python 3.12；
- [ ] `uv pip check` 没有依赖冲突；
- [ ] Python 依赖导入检查输出 `GLADE Python dependencies: OK`；
- [ ] 默认 GEM 和培养基存在；
- [ ] `chassis` 冒烟测试成功。

### 13.2 RetroPath

- [ ] 逆向规则文件存在且 SHA-256 正确；
- [ ] `docker compose ... config --quiet` 通过；
- [ ] RetroPath 镜像构建成功；
- [ ] `docker compose ... ps` 显示服务运行；
- [ ] `/health` 返回 `ready: true`；
- [ ] 服务只监听 `127.0.0.1:8765`；
- [ ] 需要预测路线验证时，MNXref `status` 通过。

### 13.3 完整设计流程

- [ ] CodonTransformer 关键模型文件存在；
- [ ] Selenzyme 地址已配置并可访问；
- [ ] Milvus 两个集合都存在且应用的 schema 检查通过；
- [ ] 需要辅助蛋白研究或文献检索时，大模型四项配置完整；
- [ ] `.env` 没有被提交到 Git；
- [ ] `inputs`、`outputs`、`model` 和 `data/retropath` 已纳入独立备份。

## 14. 常见部署问题

### 14.1 `uv sync` 使用了错误的 Python

确认版本范围并重新指定 Python 3.12：

```powershell
uv python install 3.12
uv python pin 3.12
uv sync --frozen
uv run python --version
```

输出必须是 3.12.x。

### 14.2 PowerShell 中中文显示乱码

在当前 PowerShell 会话设置 UTF-8：

```powershell
[Console]::InputEncoding = [System.Text.UTF8Encoding]::new()
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
$OutputEncoding = [System.Text.UTF8Encoding]::new()
```

然后重新运行命令。也应使用 PowerShell 7 和支持中文的终端字体。

### 14.3 找不到输入文件

确认文件确实位于 `inputs`，并且命令只传文件名：

```powershell
Get-Item inputs\deploy_smoke.json
uv run python main.py chassis -i deploy_smoke.json
```

### 14.4 Docker 无法启动或 Compose 连接失败

依次检查：

```powershell
docker version
docker info
docker compose -f compose.retropath.yml config
docker compose -f compose.retropath.yml ps
```

如果客户端存在但无法连接守护进程，通常是 Docker Desktop 尚未启动、WSL 2 后端异常，
或者当前用户没有使用 Docker 的权限。

### 14.5 RetroPath 镜像构建失败

查看失败位置后重试：

```powershell
docker compose -f compose.retropath.yml build retropath
```

常见原因包括 Docker Registry、Conda Forge、KNIME 更新站点或 Zenodo 无法访问，以及
Docker 磁盘空间不足。不要在不理解影响的情况下修改 Dockerfile 中锁定的版本。

### 14.6 RetroPath 健康检查不是 `ready: true`

```powershell
docker compose -f compose.retropath.yml logs --tail 200 retropath
Invoke-RestMethod http://127.0.0.1:8765/health | ConvertTo-Json -Depth 5
```

重点查看 `errors`。常见错误包括规则文件缺失、规则哈希不匹配、KNIME 或 RDKit 节点
缺失，以及兼容库未安装完整。规则文件修复后重新创建容器：

```powershell
docker compose -f compose.retropath.yml up -d --force-recreate retropath
```

### 14.7 RetroPath 任务因内存或超时失败

容器内存上限为 7 GiB，单任务服务超时为 1800 秒。先查看日志中的资源和终态信息。
如果 Docker Desktop 可用内存低于容器要求，应先提高 Docker 可用资源或迁移到资源更
充足的机器；不要把资源失败解释成“没有预测路线”。

### 14.8 RetroPath 验证提示缺少 MNXref

```powershell
uv run python -m src.pathway_analyze.retropath_mnxref install
uv run python -m src.pathway_analyze.retropath_mnxref status
```

如果状态提示规则哈希不同，先恢复正确的规则文件，再重建索引。

### 14.9 `protein-to-cds` 提示本地模型不完整

检查：

```powershell
Get-Item model\CodonTransformer\model.safetensors
Get-Item model\CodonTransformer\config.json
Get-Item model\CodonTransformer\tokenizer.json
```

文件缺失时重新执行第 8 节的 `hf download`。如果明确请求了 CUDA 但环境没有可用显卡，
改用：

```powershell
$inputFileName = '你的输入文件.json'
uv run python main.py protein-to-cds -i $inputFileName --device cpu
```

### 14.10 Milvus 连接成功但推荐仍失败

确认集合名称、数据库名和 schema。GLADE 不会在远端 Milvus 不可用或字段不匹配时静默
使用旧结果。若提示缺少字段或没有符合版本的数据，需要由集合管理员修复数据，不能只在
客户端修改集合名称规避检查。

### 14.11 大模型或 Selenzyme 配置失败

检查 `.env` 中相应字段是否为空、地址是否能从当前机器访问，以及 API 密钥是否有权限。
基础路线分析不依赖这些服务，可以先关闭 `--literature-search`，或暂不运行
`auxiliary-protein`，继续定位其他组件。

## 15. 更新、备份与日常运维

### 15.1 更新代码

更新前先查看本地改动：

```powershell
git status --short
```

确认可以更新后，按团队的 Git 流程拉取新版本，再重新同步锁定依赖：

```powershell
git pull --ff-only
uv sync --frozen
```

如果 `services/retropath`、`compose.retropath.yml` 或其锁定版本发生变化，重新构建并创建
服务：

```powershell
docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d --force-recreate retropath
```

更新后重新执行基础依赖检查、RetroPath 健康检查和 MNXref 状态检查。

### 15.2 需要备份的内容

至少备份：

- `.env`，使用安全的秘密管理方式保存；
- `inputs` 中的项目配置和用户上传文件；
- `outputs` 中的 `design_manifest.json` 及全部结果；
- `model/CodonTransformer` 的已审核模型快照或其精确 revision；
- `data/retropath` 中的规则与 MNXref 索引，或保留可重复安装的来源和哈希；
- RetroPath Docker 命名卷中需要长期保留的任务记录。

不要只依赖 Git 备份上述目录，因为它们大多在 `.gitignore` 中。

### 15.3 服务安全

- 保持 RetroPath 只监听 `127.0.0.1`；
- 不要把 `.env`、API 密钥或 Milvus 令牌写入文档和日志；
- 生产环境应限制项目目录、输出和缓存的访问权限；
- 使用来自指定来源且哈希正确的规则和模型；
- 不要手动编辑候选 CSV、验证 manifest 或 MNXref SQLite 来绕过一致性检查。

## 16. Linux 部署命令对照

GLADE Python 代码和 RetroPath Compose 服务可以在具备 Docker Engine 的 Linux 主机上
使用，但当前项目示例主要以 Windows 为准。Linux 上核心命令如下：

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv python install 3.12
uv sync --frozen
uv run python main.py -h

docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d retropath
curl -fsS http://127.0.0.1:8765/health | python -m json.tool
```

规则文件仍必须放在：

```text
data/retropath/rules/rr02/retrorules_rr02_rp2_flat_retro.csv
```

Linux 校验 SHA-256：

```bash
sha256sum data/retropath/rules/rr02/retrorules_rr02_rp2_flat_retro.csv
```

期望值为：

```text
e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f
```

Linux 主机同样不应把本地 RetroPath 端口改为公网监听。部署到服务器时，还需要由运维
人员根据单位规范配置用户权限、服务自启动、磁盘配额、日志轮转和秘密管理。

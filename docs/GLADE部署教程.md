# GLADE 部署教程

本文以 Windows 10/11、PowerShell 7、Python 3.12 和 `uv` 为例。所有命令均在
GLADE 项目根目录执行。部署完成后的使用方法见
[GLADE 用户使用说明](GLADE用户使用说明.md)。

docker官方下载地址：https://www.docker.com/products/docker-desktop/ 

## 1. 选择部署范围

| 功能 | 必需组件 |
|---|---|
| `chassis`、`expand`、KEGG `gap`、KEGG 路线验证 | Python 环境、GEM、培养基、KEGG 网络 |
| `gap --retropath` | 上述组件、Docker Desktop、RetroRules 文件 |
| RetroPath 路线 GEM 验证 | 上述组件、MNXref v3.0 索引 |
| 主酶、CDS、表达元件、质粒和组装 | CodonTransformer、Milvus；研究功能还需大模型和 Selenzyme |

建议先完成基础部署和冒烟测试，再安装可选组件。

基础要求：

- 64 位 Windows 10/11；
- PowerShell 7 和 Git；
- Python 3.12，不能使用 Python 3.13；
- RetroPath 需要 Docker Desktop，并使用 Linux 容器；
- RetroPath 容器上限为 7 GiB，建议主机至少有 12–16 GiB 内存；
- 基础功能需访问 Python 软件源、KEGG、Rhea 和 UniProt；其他组件还可能访问
  Docker Registry、Conda Forge、KNIME、Zenodo、MetaNetX 和 Hugging Face。

## 2. 基础部署

### 2.1 获取代码并安装依赖

进入项目根目录，确认以下文件存在：

```powershell
Get-Item main.py, pyproject.toml, uv.lock, compose.retropath.yml
```

安装 `uv`（已安装可跳过）：

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

重新打开 PowerShell，然后执行：

```powershell
uv python install 3.12
uv sync --frozen
uv run python --version
uv pip check
uv run python main.py -h
```

Python 版本应为 3.12.x。后续统一使用 `uv run python ...`，无需手动激活 `.venv`。

### 2.2 检查基础数据

默认部署使用以下文件：

```text
data/gem_models/iML1515.json
data/mediums/default_medium.json
```

检查并创建运行目录：

```powershell
Get-Item data\gem_models\iML1515.json, data\mediums\default_medium.json
New-Item -ItemType Directory -Force -Path inputs, outputs, cache, model | Out-Null
```

### 2.3 可选：配置外部服务

基础路线分析不需要 `.env`。运行辅助蛋白研究、文献检索、Selenzyme、表达元件或质粒
推荐时，在项目根目录创建 `.env`：

```dotenv
MODEL_PROVIDER=openai
AGENT_LLM_MODEL=<模型名称>
API_KEY=<API密钥>
BASE_URL=<兼容OpenAI API的服务地址>

GLADE_CONTACT_EMAIL=<联系邮箱>
SELENZYME_REST_URL=<Selenzyme REST服务地址>

MILVUS_HOST=<Milvus主机或HTTP地址>
MILVUS_PORT=19530
MILVUS_TOKEN=<可选令牌>
MILVUS_DB_NAME=<可选数据库名>
```

不要把 `.env`、API 密钥或 Milvus 令牌提交到 Git。

## 3. 可选组件

### 3.1 CodonTransformer

只有 `protein-to-cds` 需要本地模型。下载到固定目录：

```powershell
uv run hf download adibvafa/CodonTransformer --local-dir model/CodonTransformer
Get-Item model\CodonTransformer\model.safetensors
Get-Item model\CodonTransformer\config.json
Get-Item model\CodonTransformer\tokenizer.json
```

模型可在 CPU 上运行；没有正确配置 CUDA 时使用 `--device cpu` 或默认的
`--device auto`。

### 3.2 RetroPath 规则与服务

RetroPath 使用以下固定规则文件：

```text
data/retropath/rules/rr02/retrorules_rr02_rp2_flat_retro.csv
```

规则来自 [RetroRules rr02-rp2-hs](https://doi.org/10.5281/zenodo.5828017)。新克隆的
仓库通常不包含该文件。下载并安装：

```powershell
$archive = 'cache\retrorules_rr02_rp2_hs.tar.gz'
$extractDir = 'cache\retrorules_rr02_extracted'
$rulesDir = 'data\retropath\rules\rr02'

New-Item -ItemType Directory -Force -Path cache, $extractDir, $rulesDir | Out-Null
Invoke-WebRequest `
  -Uri 'https://zenodo.org/api/records/5828017/files/retrorules_rr02_rp2_hs.tar.gz/content' `
  -OutFile $archive
tar -xzf $archive -C $extractDir

$source = Get-ChildItem $extractDir -Recurse -File `
  -Filter 'retrorules_rr02_rp2_flat_retro.csv' | Select-Object -First 1
if ($null -eq $source) { throw '归档中没有找到 RetroRules 文件' }
Copy-Item $source.FullName "$rulesDir\retrorules_rr02_rp2_flat_retro.csv" -Force
```

校验 SHA-256：

```powershell
$rulesFile = 'data\retropath\rules\rr02\retrorules_rr02_rp2_flat_retro.csv'
$expected = 'e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f'
$actual = (Get-FileHash $rulesFile -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actual -ne $expected) { throw "规则文件哈希不匹配: $actual" }
```

启动 Docker Desktop 后构建并启动服务：

```powershell
docker compose -f compose.retropath.yml config --quiet
docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d retropath
Invoke-RestMethod http://127.0.0.1:8765/health | ConvertTo-Json -Depth 5
```

健康检查的 `ready` 必须为 `true`。服务仅绑定 `127.0.0.1:8765`，不要暴露到公网。

日常启动、查看日志和停止：

```powershell
docker compose -f compose.retropath.yml up -d retropath
docker compose -f compose.retropath.yml logs --tail 100 retropath
docker compose -f compose.retropath.yml down
```

不要默认使用 `down -v`，否则会删除保存任务和结果的 Docker 卷。

### 3.3 MNXref

只有验证 RetroPath 预测路线时才需要 MNXref：

```powershell
uv run python -m src.pathway_analyze.retropath_mnxref install
uv run python -m src.pathway_analyze.retropath_mnxref status
```

状态应显示 MNXref `3.0`，并且 `rr02_sha256` 与规则文件一致。规则文件变化后需要重新
构建索引。

### 3.4 Milvus

仓库只包含 Milvus 客户端，部署方需提供已填充并通过审计的集合：

| 集合 | 用途 |
|---|---|
| `expression_parts_v3` | promoter、RBS、terminator 推荐 |
| `plasmid_templates_v2` | 质粒骨架推荐 |

如果连接成功但仍提示字段或版本不匹配，应修复远端集合数据，不能仅修改集合名称绕过
检查。`model/bge-m3` 当前不参与这些推荐，无需下载。

## 4. 冒烟测试

创建 `inputs/deploy_smoke.json`：

```powershell
@'
{
  "target_name": "C00811"
}
'@ | Set-Content inputs\deploy_smoke.json -Encoding utf8NoBOM
```

`--input` 后只写 `inputs` 目录下的文件名，不要传完整路径，也不要写
`inputs/deploy_smoke.json`。

运行基础测试：

```powershell
uv run python main.py chassis -i deploy_smoke.json
uv run python main.py info -i deploy_smoke.json --chassis
uv run python main.py gap -i deploy_smoke.json -d 0
uv run python main.py info -i deploy_smoke.json --gap -d 0
```

`chassis` 成功后应生成：

```text
outputs/C00811/chassis_result/
├── producible_kegg_compounds.csv
└── analyze_chassis_metabolites_summary.csv
```

如果目标已经在底盘可生成集合中，程序会提示“目标化合物已在底盘细胞中，无需新增
合成路径”；这也是正常结果。

已部署 RetroPath 时再执行：

```powershell
uv run python main.py gap -i deploy_smoke.json --retropath -d 0
uv run python main.py info -i deploy_smoke.json --retropath -d 0
```

部署验收关注命令是否得到明确终态，不要求测试目标一定能找到候选路线。

## 5. 验收清单

- [ ] `uv run python --version` 为 Python 3.12.x；
- [ ] `uv pip check` 通过；
- [ ] `uv run python main.py -h` 可显示命令；
- [ ] GEM 和培养基文件存在；
- [ ] `chassis` 冒烟测试成功；
- [ ] 使用 RetroPath 时，规则哈希正确且 `/health` 返回 `ready: true`；
- [ ] 验证 RetroPath 路线时，MNXref `status` 通过；
- [ ] 使用 CDS 优化时，CodonTransformer 三个关键文件存在；
- [ ] 使用表达元件和质粒推荐时，两个 Milvus 集合存在；
- [ ] `.env` 未提交到 Git。

## 6. 常见问题

### Python 版本错误

```powershell
uv python install 3.12
uv python pin 3.12
uv sync --frozen
```

### 找不到输入文件

确认文件位于 `inputs`，命令只传文件名：

```powershell
Get-Item inputs\deploy_smoke.json
uv run python main.py chassis -i deploy_smoke.json
```

### PowerShell 中文乱码

```powershell
[Console]::InputEncoding = [System.Text.UTF8Encoding]::new()
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
$OutputEncoding = [System.Text.UTF8Encoding]::new()
```

### Docker 或 RetroPath 未就绪

```powershell
docker info
docker compose -f compose.retropath.yml ps
docker compose -f compose.retropath.yml logs --tail 200 retropath
Invoke-RestMethod http://127.0.0.1:8765/health | ConvertTo-Json -Depth 5
```

常见原因是 Docker Desktop 未启动、规则文件缺失或哈希错误、内存不足，以及 Docker
Registry、Conda Forge、KNIME 或 Zenodo 无法访问。

### 缺少 MNXref

```powershell
uv run python -m src.pathway_analyze.retropath_mnxref install
uv run python -m src.pathway_analyze.retropath_mnxref status
```

### CodonTransformer 不完整

重新执行第 3.1 节下载命令；没有 CUDA 时改用 `--device cpu`。

### 大模型、Selenzyme 或 Milvus 失败

检查 `.env` 是否完整、地址能否从当前机器访问以及令牌权限。基础路线分析不依赖这些
服务，可以先单独确认 `chassis` 和 KEGG `gap` 正常。

## 7. 更新与备份

更新代码：

```powershell
git status --short
git pull --ff-only
uv sync --frozen
```

如果 RetroPath 服务代码或 Compose 文件发生变化，重新构建：

```powershell
docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d --force-recreate retropath
```

需要单独备份 `.env`、`inputs`、`outputs`、`model/CodonTransformer` 和
`data/retropath`。这些内容大多不由 Git 保存。

Linux 主机使用相同的目录结构，核心命令为：

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
uv python install 3.12
uv sync --frozen
docker compose -f compose.retropath.yml up -d retropath
curl -fsS http://127.0.0.1:8765/health | python -m json.tool
```

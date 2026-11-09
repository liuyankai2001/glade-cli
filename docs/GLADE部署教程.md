# GLADE 部署教程

以 Windows 10/11、PowerShell 7、Python 3.12 和 uv 为例，命令均在项目根目录执行。日常操作见 [用户使用说明](GLADE用户使用说明.md)。

## 1. 部署范围

先完成基础部署，再按功能安装可选组件。

| 功能 | 所需组件 |
|---|---|
| 底盘分析、KEGG 路线搜索与验证 | Python 环境、GEM、培养基及网络 |
| RetroPath 路线搜索 | Docker Desktop（Linux 容器）、RetroRules、RetroPath 服务 |
| RetroPath 路线验证 | 额外安装 MNXref v3.0 |
| 从基因组构建 GEM | Bactabolize 服务、参考 GEM、参考及目标 GBFF |
| CDS 生成 | 本地 CodonTransformer 模型 |
| 表达元件与质粒推荐 | 已填充的 Milvus 集合及相应服务配置 |

准备 Git 和 PowerShell 7。需访问软件源及 KEGG、Rhea、UniProt 等数据库；可选组件还需对应下载源。Docker 安装见 [官网](https://www.docker.com/products/docker-desktop/)。使用 RetroPath 建议主机内存至少 12～16 GiB。

## 2. 基础部署

### 2.1 安装依赖

进入已有项目目录，确认包含 `main.py`、`pyproject.toml` 和 `uv.lock`。未安装 uv 时执行：

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

重新打开 PowerShell，进入项目根目录：

```powershell
cd F:\myproject\glade
uv python install 3.12
uv sync --frozen
uv run python --version
uv pip check
uv run python main.py -h
```

版本应为 3.12.x。后续使用 `uv run python ...`，无需手动激活环境。

### 2.2 准备数据

以下数据通常需随部署包单独提供：

```powershell
Get-Item data\gem_models\iML1515.json, data\mediums\default_medium.json
Get-Item data\chassis_taxonomy.json, data\codontransformer_organisms.csv
New-Item -ItemType Directory -Force -Path inputs, outputs, cache, model | Out-Null
```

前两个文件用于默认底盘分析；宿主登记表和生成类别 CSV 用于后续选酶、CDS 流程。使用新宿主时，宿主登记表所在目录需可写。

### 2.3 外部服务

基础路线分析无需 `.env`。使用大模型、Selenzyme 或 Milvus 功能时，在根目录 `.env` 中按需填写对应配置：

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

`.env` 与令牌不要提交到 Git；已有 `.env` 时合并配置，避免覆盖。

## 3. 可选组件

### 3.1 CodonTransformer

生成 CDS 前下载模型，并确认三个关键文件存在：

```powershell
uv run hf download adibvafa/CodonTransformer --local-dir model/CodonTransformer
Get-Item model\CodonTransformer\model.safetensors, model\CodonTransformer\config.json, model\CodonTransformer\tokenizer.json
```

无 CUDA 时，生成命令使用 `--device cpu`。更换宿主及生成类别见 [用户说明](GLADE用户使用说明.md#36-更换宿主)。

### 3.2 RetroPath 规则与服务

启动 Docker Desktop，选择 Linux 容器。下载 [RetroRules rr02-rp2-hs](https://doi.org/10.5281/zenodo.5828017) 并安装规则：

```powershell
$archive = 'cache\retrorules_rr02_rp2_hs.tar.gz'
$extractDir = 'cache\retrorules_rr02_extracted'
$rulesDir = 'data\retropath\rules\rr02'
New-Item -ItemType Directory -Force -Path cache, $extractDir, $rulesDir | Out-Null
Invoke-WebRequest -Uri 'https://zenodo.org/api/records/5828017/files/retrorules_rr02_rp2_hs.tar.gz/content' -OutFile $archive
tar -xzf $archive -C $extractDir
$source = Get-ChildItem $extractDir -Recurse -File -Filter 'retrorules_rr02_rp2_flat_retro.csv' | Select-Object -First 1
if ($null -eq $source) { throw '归档中没有找到 RetroRules 文件' }
Copy-Item $source.FullName "$rulesDir\retrorules_rr02_rp2_flat_retro.csv" -Force
$actual = (Get-FileHash "$rulesDir\retrorules_rr02_rp2_flat_retro.csv" -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actual -ne 'e24eb97d3172195d03abed6e7da07a4cfd53965553853d126aaa8a93b4bc552f') { throw '规则文件哈希不匹配' }
```

构建并启动：

```powershell
docker info
docker compose -f compose.retropath.yml config --quiet
docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d retropath
Invoke-RestMethod http://127.0.0.1:8765/health | ConvertTo-Json -Depth 5
```

等待 `ready: true` 后使用。默认容器硬上限 8 GiB、任务软上限 7.5 GiB；可在 `.env` 调整：

```dotenv
RETROPATH_CONTAINER_MEMORY_LIMIT=8g
RETROPATH_MEMORY_LIMIT_BYTES=8053063680
```

软上限须小于硬上限，硬上限须低于 Docker/WSL 可用内存。修改后执行 `docker compose -f compose.retropath.yml up -d --force-recreate retropath`，检查健康响应中的 `memory_limit_bytes`。

服务仅绑定本机 `127.0.0.1:8765`。详细运行与排错见 [RetroPath 本地服务说明](RetroPath本地服务使用说明.md)。

### 3.3 MNXref

仅验证 RetroPath 路线时需要：

```powershell
uv run python -m src.pathway_analyze.retropath_mnxref install
uv run python -m src.pathway_analyze.retropath_mnxref status
```

确认版本为 `3.0`、索引状态正常。规则文件变化后需重新构建索引。

### 3.4 Milvus

仓库提供客户端，部署方需准备已填充并通过审计的集合：

| 集合 | 用途 |
|---|---|
| `expression_parts_v3` | 启动子、RBS、终止子推荐 |
| `plasmid_templates_v2` | 质粒骨架推荐 |

连接配置见第 2.3 节。字段或版本不匹配时修复远端数据；`model/bge-m3` 无需下载。

### 3.5 Bactabolize 模型构建服务

已有 GEM 可跳过。需要从注释基因组构建模型时，先启动独立服务：

```powershell
docker compose -f compose.bactabolize.yml config --quiet
docker compose -f compose.bactabolize.yml build bactabolize
docker compose -f compose.bactabolize.yml up -d bactabolize
Invoke-RestMethod http://127.0.0.1:8766/health | ConvertTo-Json -Depth 5
```

等待 `ready: true`。默认 2 CPU、2 GiB 内存，可在 `.env` 设置 `BACTABOLIZE_CONTAINER_CPUS`、`BACTABOLIZE_CONTAINER_MEMORY_LIMIT`，再用 `up -d --force-recreate bactabolize` 应用。容器使用独立锁定环境，无需调整宿主 Python 版本。

准备参考 GEM（COBRA JSON）、与参考模型基因 ID 对应的参考 GBFF、已注释的目标 GBFF；裸 FASTA 不适用。创建 `inputs/gem_build_demo.json`，替换示例路径：

```json
{
  "target_name": "C00811",
  "model_file": "outputs/gem_models/my_strain.json",
  "bactabolize": {
    "reference_model_file": "data/gem_models/reference.json",
    "reference_genome_file": "data/genomes/reference.gbff",
    "target_genome_file": "data/genomes/my_strain.gbff"
  }
}
```

```powershell
uv run python main.py gem-build -i gem_build_demo.json
uv run python main.py chassis -i gem_build_demo.json
```

`model_file` 为必填输出路径；三个嵌套文件为必填输入，路径相对项目根目录或使用绝对路径。默认相似度 80%、覆盖度 25%、M9、有氧；需调整时使用 `min_pident`、`min_coverage`、`media_type`、`atmosphere_type`。参考模型无法唯一推导生物量目标时，补充 `biomass_reaction_id`。

每次 `gem-build` 提交新任务，仅保存模型 JSON；检查终端任务状态和 FBA 复验提示，未通过生长验证的模型只能作为草稿。后续选酶需另设顶层 `genome_file`，不会自动使用嵌套的目标 GBFF。自定义宿主配置见 [用户说明](GLADE用户使用说明.md#36-更换宿主)。

## 4. 部署验收

创建 UTF-8 文件 `inputs/deploy_smoke.json`：

```json
{
  "target_name": "C00811"
}
```

**`-i` / `--input` 只跟文件名**，不写 `inputs/` 或完整路径。

```powershell
uv run python main.py chassis -i deploy_smoke.json
uv run python main.py info -i deploy_smoke.json --chassis
uv run python main.py gap -i deploy_smoke.json -d 0
uv run python main.py info -i deploy_smoke.json --gap -d 0
```

确认 `outputs/C00811/chassis_result/` 下有 `producible_kegg_compounds.csv` 和 `analyze_chassis_metabolites_summary.csv`。目标已可供给或搜索无路线都可能是正常结果，验收不要求一定找到路线。

已安装 RetroPath 时再运行：

```powershell
uv run python main.py gap -i deploy_smoke.json --retropath -d 0
uv run python main.py info -i deploy_smoke.json --retropath -d 0
```

基础验收需满足 Python 3.12.x、`uv pip check` 通过、CLI 帮助可用和底盘分析成功。可选组件分别确认模型文件、健康状态、MNXref 索引或 Milvus 集合；Bactabolize 还需实际提交构建并检查复验提示。

## 5. 常见问题

| 问题 | 处理 |
|---|---|
| Python 版本错误 | 执行 `uv python install 3.12`、`uv python pin 3.12`、`uv sync --frozen` |
| 找不到输入 | 文件放在 `inputs`，`-i` 只传文件名 |
| Docker 无法连接 | 启动 Docker Desktop，执行 `docker info` 检查 Linux 容器环境 |
| 服务未就绪或任务失败 | 检查容器状态、日志、资源和下载源；Bactabolize 还需核对模型与 GBFF 基因 ID |
| 缺少 MNXref 或模型文件 | 重做第 3.3 或 3.1 节；无 CUDA 时使用 `--device cpu` |
| 大模型、Selenzyme、Milvus 失败 | 检查 `.env` 中的地址、令牌权限及对应服务状态 |

PowerShell 中文乱码时执行：

```powershell
[Console]::InputEncoding = [System.Text.UTF8Encoding]::new()
[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
$OutputEncoding = [System.Text.UTF8Encoding]::new()
```

## 6. 维护与备份

服务检查命令如下；Bactabolize 将文件名和服务名分别换成 `compose.bactabolize.yml`、`bactabolize`：

```powershell
docker compose -f compose.retropath.yml ps
docker compose -f compose.retropath.yml logs --tail 100 retropath
docker stats --no-stream
```

更新代码前检查本地改动，再更新依赖：

```powershell
git status --short
git pull --ff-only
uv sync --frozen
```

服务代码或构建配置更新后，重新 `build` 并 `up -d --force-recreate`。停止服务用对应 Compose 的 `down`；`down -v` 会删除任务和结果卷。

单独备份 `.env`、`inputs/`、`outputs/`、`model/CodonTransformer/` 和所用 `data/` 文件（包括 GEM、培养基、宿主登记表、生成类别 CSV、RetroPath 数据）；需要服务历史时另备份 Docker 作业卷。这些内容大多不随 Git 保存。

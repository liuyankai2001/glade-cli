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
| 可选 `gem-build`（从注释基因组构建 GEM） | Docker Desktop、Bactabolize 服务、参考 GEM JSON、参考与目标 GBFF |
| 主酶、CDS、表达元件、质粒和组装 | CodonTransformer、Milvus；研究功能还需大模型和 Selenzyme |

建议先完成基础部署和冒烟测试，再安装可选组件。

基础要求：

- 64 位 Windows 10/11；
- PowerShell 7 和 Git；
- Python 3.12，不能使用 Python 3.13；
- RetroPath 需要 Docker Desktop，并使用 Linux 容器；
- RetroPath 容器默认上限为 7 GiB，可通过项目根目录 `.env` 调整；建议主机至少有
  12–16 GiB 内存；
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

基础路线分析不需要 `.env`。运行 Selenzyme、表达元件或质粒
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

RetroPath 内存限制可以写入项目根目录已有的 `.env`。模板见 `.env.example`：

```dotenv
# Docker 容器硬上限
RETROPATH_CONTAINER_MEMORY_LIMIT=8g

# RetroPath 工作集软上限，单位为字节（这里是 7.5 GiB）
RETROPATH_MEMORY_LIMIT_BYTES=8053063680
```

软上限必须低于容器硬上限，并且容器硬上限不能超过 Docker/WSL 实际可用内存。修改后
需要重建容器才能生效。可以先检查 Docker 总内存：

```powershell
docker info --format '{{.MemTotal}}'
```

启动 Docker Desktop 后构建并启动服务：

正式构建使用项目根目录作为上下文，仅发送 `services/retropath` 和指定 RR02 规则文件；Dockerfile 专用忽略规则排除 `.env`、实验输出及其他项目内容。镜像从 micromamba 基础镜像安装锁定运行栈，构建时编译源码补丁工具和 agent；不依赖旧 RetroPath 镜像或 docs 中的预编译 JAR。

新镜像标签为 `glade-retropath:3.9.1-r20260212-memoryfix2`。规则与修补内置镜像，运行时只挂载作业数据卷；无需旧的 8g/memoryfix Compose 覆盖配置。默认 Java 堆与健康检查均为 1024 MiB、CPU 限额为 1，内存加交换空间上限为 13 GiB。

```powershell
docker compose -f compose.retropath.yml config --quiet
docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d retropath
Invoke-RestMethod http://127.0.0.1:8765/health | ConvertTo-Json -Depth 5
```

修改 `.env` 中的内存配置后执行：

```powershell
docker compose -f compose.retropath.yml up -d --force-recreate retropath
```

健康检查中的 `memory_limit_bytes` 应等于 `.env` 中配置的软上限。Docker Desktop 使用
WSL 2 时，WSL 总内存仍由用户目录中的 `.wslconfig` 单独控制。

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

### 3.5 Bactabolize 模型构建服务（可选）

已有 GEM 时可直接运行 `chassis`，无需安装本组件。需要从目标菌株基因组构建
GEM 时，先部署 Bactabolize，再调用 `gem-build`。服务独立于 RetroPath，默认地址为
`http://127.0.0.1:8766`，两个服务可以分别启动。

本组件要求 Docker Desktop 使用 Linux 容器。正式镜像固定 Bactabolize 1.0.6、
Python 3.9.0、COBRApy 0.21、Biopython 1.79、BLAST 2.16 和 Prodigal 2.6.3；
科学运行栈由 `services/bactabolize/environment-linux-64.lock` 固定，HTTP 层单独锁定。
宿主项目仍使用原有 Python 3.12 和 `uv.lock`。构建上下文仅为
`services/bactabolize`，不包含项目 `.env`、数据或实验目录，也不依赖 `tests`。
该历史科学环境的部分上游依赖声明与实际锁定版本不一致（例如 COBRApy 的
pandas、Pydantic 声明），全局 `pip check` 会报告既有冲突；不要据此升级科学运行栈。
镜像构建单独校验核心版本和新增 HTTP 包的依赖，实际模型生成仍需结合验证报告验收。

首次部署（在项目根目录执行）：

```powershell
docker info
docker compose -f compose.bactabolize.yml config --quiet
docker compose -f compose.bactabolize.yml build bactabolize
docker compose -f compose.bactabolize.yml up -d bactabolize
Invoke-RestMethod http://127.0.0.1:8766/health | ConvertTo-Json -Depth 5
```

等待 `/health` 的 `ready` 为 `true`，检查 `versions`。首次冷构建需要从 Conda
和 PyPI 下载锁定的软件包，耗时取决于网络；已有镜像后的容器启动无需重复安装。
构建和启动属于部署步骤，日常每次 `gem-build` 只提交任务，不启动 Docker 服务。

服务以非 root 用户运行，单工作进程串行执行构建，队列总容量为 8（含运行中任务），
单任务最长 1800 秒，上传请求上限为 64 MiB。作业和 SQLite 状态保存在
`glade-bactabolize_bactabolize-jobs` 命名卷；重启后排队任务继续处理，已中断的
运行中任务标为失败，不会被当作成功结果。

默认限制为 2 CPU、2 GiB 内存。可在项目根目录已有 `.env` 中加入：

```dotenv
BACTABOLIZE_CONTAINER_CPUS=2.0
BACTABOLIZE_CONTAINER_MEMORY_LIMIT=2g
```

调整后执行 `docker compose -f compose.bactabolize.yml up -d --force-recreate bactabolize`。
资源限制应小于 Docker/WSL 实际可用资源。服务内部关闭宽限为 30 秒，Compose 等待
40 秒后才强制终止，以便清理子进程。端口只绑定本机；不要直接暴露到公网。

准备三个文件：参考 GEM（COBRA JSON）、与参考模型基因 ID 对应的参考 GBFF、
已注释的目标 GBFF。此流程固定使用已有注释（`no_reannotation=true`），不接受裸
FASTA 替代 GBFF。下例路径是示例，须换成实际文件；数据文件不会随镜像提供。
创建 `inputs/gem_build_demo.json`：

```json
{
  "target_name": "C00811",
  "model_file": "outputs/gem_models/my_strain.json",
  "bactabolize": {
    "reference_model_file": "data/gem_models/reference.json",
    "reference_genome_file": "data/genomes/reference.gbff",
    "target_genome_file": "data/genomes/my_strain.gbff",
    "min_pident": 80,
    "min_coverage": 25,
    "media_type": "m9",
    "atmosphere_type": "aerobic",
    "biomass_reaction_id": null
  }
}
```

JSON 内文件路径相对于项目根目录，也可使用绝对路径。`gem-build` 必须明确提供
`model_file` 作为输出，且不要求 `target_name`；例中加入目标是为了接着运行
`chassis`。三个嵌套输入文件必填，其余 Bactabolize 参数均可省略，默认值如上。
相似度与覆盖度取值为 `(0, 100]`；培养基支持 `m9`/`lb`，气氛支持
`aerobic`/`anaerobic`。生物量 ID 为 `null` 时从参考模型唯一正向目标反应推导，
无法唯一确定时应显式指定。`medium_path` 可在顶层指定 `gem-build` 的宿主 FBA
复验培养基，默认 `data/mediums/default_medium.json`，与服务的 `media_type`
分别使用；该选项仅用于构建复验，后续分析命令继续使用项目已有培养基配置。

```powershell
uv run python main.py gem-build -i gem_build_demo.json
uv run python main.py chassis -i gem_build_demo.json
```

`--input`/`-i` 只跟 `inputs` 下文件名。普通分析命令的 `model_file` 省略、`null`
或空白时仍使用默认 `data/gem_models/iML1515.json`；`gem-build` 要求显式输出路径。
顶层 `genome_file` 用于已有主酶流程的基因组 EC 过滤，与这里的目标 GBFF 输入独立，
不会自动由 `bactabolize.target_genome_file` 替换。本次不扩展蛋白选型和 CDS 的
Nissle 1917 宿主预设。

自定义模型的分析输出按模型内容隔离，例如
`outputs/C00811/model_<SHA256前16位>/`。修改模型内容后会使用新的结果目录，
需要重新运行 `chassis` 再继续路线分析；默认 iML1515 保持原有输出目录。

每次执行都会提交新的重建任务，已有 JSON 不会触发缓存复用。构建完成后只保存
`model_file` 指定的原始模型 JSON；输出目录中不新增 XML、`model_build`、日志或
验证报告，也不向模型添加缓存字段。原始 Bactabolize 注释和救援基因别名保留。

客户端在系统临时目录中完成输入/结果哈希核验、JSON 结构检查、COBRApy 读取和
FBA 复验，清理临时数据后原子发布模型；终端显示任务 ID、模型路径和必要提示。
并发锁存放在系统临时目录；输出请使用独立子目录，不能把包含锁目录的系统临时
根目录直接作为输出目录，否则会在提交任务前报 `output_temp_conflict`。
生成可读模型与验证能否生长是不同结果：生长失败时模型仍可发布为未验证草稿，
终端会提示；它不能证明模型具有生物学准确性。任务失败、超时、取消或结果
校验失败不会覆盖已有正式模型。本地不保留任务恢复记录，排查使用终端任务 ID
及 Docker 服务内的日志和结果。

旧 `--force` 参数继续接受，与普通运行行为相同。之前生成的 XML 和 `model_build`
不会被自动删除，但后续运行不再读写这些产物。服务地址可用嵌套
`bactabolize.service_url` 或环境变量 `BACTABOLIZE_SERVICE_URL` 覆盖；默认使用本机 8766。

日常维护与故障排查：

```powershell
docker compose -f compose.bactabolize.yml ps
docker compose -f compose.bactabolize.yml logs --tail 200 bactabolize
Invoke-RestMethod http://127.0.0.1:8766/health | ConvertTo-Json -Depth 5
docker stats --no-stream
docker compose -f compose.bactabolize.yml down
```

`down` 保留作业卷；不要默认使用 `down -v`，它会删除任务、上传文件和结果。
健康检查失败先检查软件版本、卷写权限和服务日志；无法连接时检查 Docker、容器状态
和端口；构建失败检查 Conda/PyPI 网络；作业失败根据终端的任务 ID 查询服务任务
接口及容器任务卷内日志，并检查模型/GBFF 基因 ID 对应关系。出现 OOM、超时或
队列已满时按诊断调整资源或等待空闲后重试。服务代码或锁文件更新后重新 `build`
并 `up -d --force-recreate`；备份模型 JSON，需要服务端历史时另行备份作业卷。

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
- [ ] 使用 `gem-build` 时，Bactabolize `/health` 就绪，构建产物与验证报告已检查；
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

## 8. 可选用法速查

以下命令使用主文档中的 `demo01.json`。使用 uv 时，将 `python main.py` 换成
`uv run python main.py`。常规完整流程见 [用户使用说明](GLADE用户使用说明.md)。

### 路线验证

指定多条路线；省略 `-s` 时验证当前深度的全部路线：

```powershell
python main.py validate -i demo01.json -s 1 2 3 -m per -c strict -d 0
```

| 参数 | 用途 |
|---|---|
| `-m per` | 独立验证每条路线，默认模式 |
| `-m pooled` / `-m both` | KEGG 路线联合验证 / 独立加联合验证 |
| `-c strict` / `-c relaxed` | 严格辅因子模式 / 放宽辅因子以诊断阻断 |

包含 RetroPath 路线时使用 `-m per`，并先完成第 3.3 节的 MNXref 安装。
验证、查看和写入使用与搜索相同的 `-d`。未验证或验证失败仍可选择，但需要人工复核。

查看路线更多信息：

```powershell
python main.py info -i demo01.json --solution 1 --step 2 -d 0
python main.py info -i demo01.json --solution 1 --all -d 0
python main.py info -i demo01.json --solution 1 --all --verbose -d 0
```

`--all` 展示培养基有机底物到目标的完整路线，需要保留模型、培养基和底盘结果，
不能与 `--step` 同用。`--verbose` 进一步显示反应和通量详情。

### 更多可选参数

| 目的 | 命令示例 |
|---|---|
| 每步保留 10 个主酶候选 | `python main.py main-enzyme -i demo01.json --top-n 10` |
| 最多生成 10 个主酶组合 | `python main.py main-enzyme-sets -i demo01.json --max-sets 10` |
| 查看第 1 步的主酶候选 2 | `python main.py info -i demo01.json --main-enzyme-candidate 2 --step 1` |
| 使用 CPU 生成 CDS | `python main.py protein-to-cds -i demo01.json --device cpu` |
| 请求 24 个元件方案 | `python main.py expression --design --parts -i demo01.json --n-designs 24` |
| 请求 10 个质粒候选 | `python main.py plasmid --recommend -i demo01.json --n-candidates 10` |
| 优先卡那霉素抗性 | `python main.py plasmid --recommend -i demo01.json --preferred-resistance kanamycin` |
| 排除指定抗性 | `python main.py plasmid --recommend -i demo01.json --exclude-resistance ampicillin tetracycline` |
| 指定双酶切 / Gibson 组装 | `python main.py assembly --plan -i demo01.json --method restriction` / `--method gibson` |

元件方案数量范围 3～96、默认 12；质粒候选数量范围 1～20、默认 5。

主酶选型先使用 UniProt、Rhea、KO 等常规数据库，仍缺候选时使用 Selenzyme 兜底。
自动搜索论文及抽取酶活证据的功能已移除；数据库自带的论文引用、实验注释和
交叉引用继续保留。旧 `--literature-search` 参数不再接受，也不再生成文献检索审计文件。
升级前的酶候选、组合及 manifest 主酶选择结果需要重新生成：运行 `main-enzyme`、
`main-enzyme-sets`，再通过 `write --main-enzyme-set N` 写入组合，并重新构建下游设计。
原有文件不会自动删除；此次升级无需重新运行 `gem-build`、`chassis` 或路线搜索。
质粒优先级支持 `--priority stability`、`balanced`、`expression`，默认 `stability`。
指定组装方法后，任一方案不可行都会使计划不完整，需要处理后才能接受。

CDS 生成阶段还可附加 motif 统计：

```powershell
python main.py protein-to-cds -i demo01.json --forbidden-motif GAATTC --forbidden-motif GGATCC
```

这个参数只统计，实际位点消除使用 `optimize --enzyme`。
重新生成 CDS 会覆盖 optimized，已有编辑需要重新执行。

# 表达盒分组与选择

完成密码子优化后，先生成系统推荐的表达盒分组方案：

```powershell
uv run python main.py expression --design --box -i demo04.json
```

命令会列出可用的 `design_id`，并将完整候选结果写入：

```text
outputs/<目标化合物>/expression_box/expression_box_designs.json
```

确认方案后，将对应编号写入 `design_manifest.json`。例如选择方案 1：

```powershell
uv run python main.py write -i demo04.json --expression-box 1
```

写入结果位于 manifest 的 `expression_box_selection`。相同方案重复写入不会增加
manifest revision；如果 CDS 或更上游的路线、蛋白选择发生变化，旧表达盒选择会被
自动清除，需要重新生成并选择方案。

该命令只确认蛋白如何分组，不会选择 promoter、RBS、terminator，也不会组装
GenBank 文件。

## 推荐表达元件

表达盒分组写入 manifest 后，可以从远端 Milvus 推荐表达元件：

```powershell
uv run python main.py expression --design --parts -i demo04.json
```

默认生成 12 个方案。也可以根据实验容量指定 3 到 96 个方案，例如：

```powershell
uv run python main.py expression --design --parts -i demo04.json --n-designs 24
```

命令读取 `.env` 中的 `MILVUS_HOST` 和 `MILVUS_PORT`，从
`expression_parts_v3` 获取适用于 MG1655 的 promoter、RBS 和 terminator。
它使用结构化字段过滤，不需要 BGE-M3、Torch 或本地 embedding 模型。

系统不再强制平均分配低负担、均衡和高表达方案。所有通过安全检查的组合使用统一的
稳定表达评分排序，分数最高的前 N 个入选；低、中、高仅作为表达水平标签。

稳定表达评分为 0 到 100 分，综合元件证据与宿主匹配、promoter 稳健性、各 CDS
的 OSTIR 翻译稳健性、terminator 可靠性和重复元件风险。它是可解释的设计启发式，
不是实验成功概率，也不能证明酶活性或目标化合物产量。

每个方案还会独立计算 `expression_burden.v1` 数值负担。系统把 promoter 活性百分位、
RBS 活性百分位、OSTIR 在同一基因上下文候选中的经验百分位、CDS 长度和非预期翻译
起始位点组合为 0–100 分：低于 35 为 `low`，35 至 65 为 `moderate`，65 及以上为
`high`。这个分数用于后续动态匹配质粒拷贝类型，不参与当前稳定表达成功分排序。
完整候选保存逐基因计算过程，因此可以追溯每一项负担来源。

每个 RBS 都会结合对应优化 CDS 的上下文，用 OSTIR 重新计算翻译起始速率；
随后使用 DNA Chisel 检查完整表达盒的 GC、同聚物和禁用酶切位点。结果写入：

```text
outputs/<目标化合物>/expression_box/expression_parts_designs.json
```

终端只显示方案数量、70 分入选门槛、分数范围、主推荐 `design_id` 及其负担摘要；完整元件组合
和分项评分仅写入上述 JSON 文件，避免大量序列细节占满命令行。

只有评分不低于 70 分的唯一安全组合才会输出。如果合格组合少于请求数量，文件保存
全部合格方案并标记 `partial`；没有方案达到门槛时标记 `failed`。系统不会使用低分
或重复方案凑数。

推荐命令本身不会修改 manifest。确认要保留的构建方案后，可以一次写入单个、多个或
闭区间编号：

```powershell
uv run python main.py write -i demo04.json --expression-parts 1
uv run python main.py write -i demo04.json --expression-parts 1 3 5
uv run python main.py write -i demo04.json --expression-parts 1:12
uv run python main.py write -i demo04.json --expression-parts 1:4 7 9:12
```

`start:end` 包含两端；重复编号自动去重，最终按候选 `rank` 排序，排名最高的已选
方案成为 `primary_design_id`。manifest 的 `parts_selection` 保存方案 ID、评分、紧凑的
数值负担摘要和内容指纹，不复制完整元件组合；候选文件发生变化后，下游必须拒绝使用旧引用。

`expression_parts_designs.v3` 和 `parts_selection.v1` 不包含完整负担数据，不能用于动态
质粒评分。升级后需要重新运行表达元件推荐和写入命令。

写入选择的同时，系统会为每个已选方案自动生成一份完整串联 GenBank：同一方案内
的表达盒按 `cassette_index` 排序后直接首尾连接，不插入 linker；每个表达盒内部按
`promoter -> RBS -> CDS -> ... -> terminator` 复原。选择 `1:12` 会生成：

```text
outputs/<目标化合物>/expression_constructs/design_001.gb
...
outputs/<目标化合物>/expression_constructs/design_012.gb
```

每个 GenBank 包含 source、表达盒范围、promoter、RBS、CDS 和 terminator feature，
manifest 的 `assembled_expression_constructs` 保存文件哈希、完整序列哈希、盒坐标和
整段安全审计。系统会先核对每个表达盒的长度与哈希，再对串联后的完整序列重新检查
GC、同聚物和禁用酶切位点，因此盒与盒连接处新产生的风险也会阻止整批写入。

文件采用整批暂存和替换：任一方案失败时不修改 manifest，也不留下部分 GenBank。
同样的选择重复执行时会复用完整文件并保持 manifest revision 不变；文件丢失或损坏
时会自动重建整批文件，但如果记录内容未变化也不会增加 revision。更换选择后，固定
的 `expression_constructs` 目录只保留本次选择对应的文件。

远端 Milvus 不可达、collection schema 不匹配或缺少任一类可用元件时，推荐命令会
直接失败，不使用本地或过期缓存静默回退。

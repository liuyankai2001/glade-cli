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

每个 RBS 都会结合对应优化 CDS 的上下文，用 OSTIR 重新计算翻译起始速率；
随后使用 DNA Chisel 检查完整表达盒的 GC、同聚物和禁用酶切位点。结果写入：

```text
outputs/<目标化合物>/expression_box/expression_parts_designs.json
```

终端只显示方案数量、70 分入选门槛、分数范围和主推荐 `design_id`；完整元件组合
和分项评分仅写入上述 JSON 文件，避免大量序列细节占满命令行。

只有评分不低于 70 分的唯一安全组合才会输出。如果合格组合少于请求数量，文件保存
全部合格方案并标记 `partial`；没有方案达到门槛时标记 `failed`。系统不会使用低分
或重复方案凑数。

这个阶段不会修改 manifest，也不会生成最终 GenBank。远端 Milvus 不可达、
collection schema 不匹配或缺少任一类可用元件时，命令会直接失败，不使用本地或
过期缓存静默回退。

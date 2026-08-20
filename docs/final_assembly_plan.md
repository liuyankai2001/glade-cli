# 最终组装计划

`assembly --plan` 为同一个已选质粒骨架和每个完整表达构建分别生成一条组装计划。它只决定组装方法、插入位置、骨架线性化方式和限制酶，不生成最终质粒 GenBank/FASTA。

## 自动推荐

```powershell
uv run python main.py assembly --plan -i demo04.json
```

系统对每个 `expression_constructs/design_NNN.gb` 独立比较 Gibson 和双酶切，并且只保留该 design 的最高分可行方案。不同 design 可以采用不同 method、酶和插入位置。

终端仅显示：design ID、method、分数、插入/替换坐标、骨架线性化方式，以及限制酶名称或 `none (PCR linearization)`。完整参数写入：

```text
outputs/<目标化合物>/final_assemble_plan/assembly_plan_recommendations.json
```

## 用户统一指定 method

```powershell
uv run python main.py assembly --plan -i demo04.json --method restriction
uv run python main.py assembly --plan -i demo04.json --method gibson
```

用户只统一指定 method；系统仍然为每个 design 独立确定插入位置、酶对或 Gibson 线性化方式。指定 method 对任一 design 不可行时，结果标记为 `partial`，列出失败 design，并禁止写入 manifest；系统不会偷偷改用另一种 method。

## 酶和线性化字段

双酶切计划记录：

- 左右限制酶名称；
- 识别序列；
- 识别位点和切割坐标；
- 替换范围；
- `restriction_site_retention=retain`。

Gibson 会比较两种骨架线性化方式：

- `restriction`：记录一个开环酶，或 cargo replacement 时记录左右两个酶；
- `pcr`：`restriction_enzymes=[]`，摘要固定为 `none (PCR linearization)`。

第一版使用 30 bp 同源臂；同源臂必须在环状骨架中唯一、GC 为 30–70%，且最大同聚物不超过 7 bp。限制酶来自常用酶目录，切位点由 Biopython `Bio.Restriction` 验证。酶切方案要求两种酶在骨架各切一次、在当前 insert 中均不切、位点处于同一个 audited insertion region，并且不与可定位的 protected features 重叠。SEVA cargo replacement 固定使用 PacI/SpeI。

## 接受计划

完整的 12 项计划生成后运行：

```powershell
uv run python main.py write -i demo04.json --assembly-plan
```

该命令将整套计划写入 manifest 的 `final_assembly_plan.v2`。相同计划重复写入不增加 revision；计划变化会清理旧 `final_assembly` 和 `final_design_report`。

`final_assembly_plan` 包含 12 个 `design_plans`，每项保存：

- `parts_design_id` 和 insert 文件哈希；
- backbone 文件哈希和 assembly policy；
- `assembly_method`；
- `target` 插入/替换坐标；
- `backbone_linearization` 及所用酶；
- restriction 或 Gibson 参数；
- 预计最终长度、评分、理由、警告和稳定指纹。

## 与执行工具的边界

规划阶段不会拼接最终序列、重排最终 features 或输出 `.gb/.fasta`。后续接入执行工具后，由执行命令遍历 manifest 中 12 个 `design_plans`，完成精确序列组装和输出校验。

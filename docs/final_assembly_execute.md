# 最终理论组装执行

`assembly --execute` 读取已经接受的 `final_assembly_plan.v2`，为每个表达元件方案生成一个完整的理论质粒。它不会重新推荐组装方法、插入位置或限制酶。

## 使用顺序

```powershell
# 1. 由 final_assemble_plan 生成计划
uv run python main.py assembly --plan -i demo04.json

# 2. 用户接受完整计划
uv run python main.py write -i demo04.json --assembly-plan

# 3. 由 final_assemble_execute 执行计划
uv run python main.py assembly --execute -i demo04.json
```

不提供单独的 `assembly --validate`。执行命令内部会核对 manifest schema、计划指纹、骨架和 insert 文件哈希、序列哈希、限制酶位点、Gibson 同源臂以及输出序列一致性。

## 输出

默认写入：

```text
outputs/<目标化合物>/final_assembly/
├── design_001_final.gb
├── design_001_final.fasta
├── design_001_assembly.json
├── ...
├── design_012_final.gb
├── design_012_final.fasta
├── design_012_assembly.json
├── run_summary.json
└── final_design_report_zh.md
```

GenBank 会保留并重新定位骨架和完整表达构建中的 features；每个 JSON 记录来源计划、method、酶、插入位置、文件哈希和计算验证结果。

中文报告在执行时由 Python 固定模板自动生成，不调用 LLM，也不读取模型配置。报告中的通路、KEGG Orthology、宿主、主酶风险、优化 CDS、表达方案、质粒和组装结果都直接映射自当前 manifest 与输出文件，`generated_by` 固定为 `system_template`。

报告固定区分：

- KO 表示 KEGG Orthology，不表示 gene knockout；
- 缺少宿主敲除分析时明确写为“暂无敲除建议”；
- 反应记录解析完成不等于主酶反应适配已经验证；
- 已生成的优化 CDS 会显示文件路径，不会再次建议从零做密码子优化；
- 复制子家族和拷贝数等级分别显示。

## complete、partial 和 failed

- `complete`：全部计划成功生成；
- `partial`：部分 design 失败，成功 design 的文件仍会保留；
- `failed`：没有 design 成功，但仍会保存运行摘要和失败报告。

每次重新执行都会在临时 staging 目录中重新生成整套结果，然后替换旧的 `final_assembly` 目录，不会混用上一次的成功文件。目录安装和 manifest 更新是一个可回滚事务。

执行结果写入 manifest：

- `final_assembly.v2`：成功构建、失败条目、文件和 bundle 指纹；
- `final_design_report.v2`：中文报告文件、生成方式和来源 bundle 指纹。

所有输出均为理论计算设计，不代表已经完成湿实验组装、转化、表达或目标化合物合成。

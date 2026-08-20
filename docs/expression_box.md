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

# 从 manifest 生成优化 CDS

`protein-to-cds` 直接读取当前项目的 `design_manifest.json`，无需用户再次输入蛋白 accession。

## 使用

先同步项目环境：

```powershell
uv sync
```

执行：

```powershell
uv run python main.py protein-to-cds -i config.json
```

默认自动选择 CUDA 或 CPU。也可以明确指定设备，并重复添加额外禁用 motif：

```powershell
uv run python main.py protein-to-cds -i config.json --device cpu `
  --forbidden-motif GAATTC --forbidden-motif GGATCC
```

本流程不需要 `torchvision` 或 `torchaudio`。

## 输入规则

- 主酶来自 `main_enzyme_selection.proteins`。
- 若存在 `auxiliary_protein_selection`，只有 `auxiliary_proteins_to_introduce` 中的蛋白会被加入，并从 `main_enzymes[*].confirmed_auxiliary_proteins` 读取详情。
- `auxiliary_protein_selection.can_advance` 必须为 `true`；该区段不存在时只处理主酶并写入警告。
- 目前仅支持 `ecoli_mg1655`，映射到 CodonTransformer organism ID 52。
- 主酶写入 manifest 时不再记录氨基酸序列哈希；运行本流程时直接读取或下载 FASTA。

## 下载、缓存与输出

输出目录为：

```text
outputs/<target>/protein_to_cds/
├── protein_sequences/<accession>.fasta
├── raw_cds/<accession>.raw.fasta
├── optimized_cds/<accession>.optimized.fasta
├── reports/<accession>.optimization.json
└── run_summary.json
```

如果本地蛋白 FASTA 是合法、单条且 accession 一致的记录，流程会直接复用，不访问 UniProt。无效或 accession 冲突的文件不会被覆盖。通过约束门禁且输入指纹一致的优化结果也会直接复用。

每个最终 CDS 都必须保持蛋白翻译完全一致，并通过起止密码子、全局及局部 GC、CAI、稀有密码子簇、双链禁用 motif 和同聚物检查。

## manifest 与退出状态

运行结果会整体写入 `cds_selection`：

- `complete`：全部成功，CLI 退出码为 0。
- `partial`：部分成功；成功产物和失败原因都会写入 manifest，CLI 退出码为 2。
- `failed`：全部失败；所有失败原因都会写入 manifest，CLI 退出码为 2。

处理期间若 manifest revision 发生变化，本次结果不会覆盖新版本。重新选择路线、主酶或辅助蛋白研究结果时，旧 `cds_selection` 会被自动清除。

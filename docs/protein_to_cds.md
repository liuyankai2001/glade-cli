# 从 manifest 生成或接收 CDS

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

## 手动添加辅助序列

把 FASTA、FAA 或纯文本文件放入 `inputs`，然后明确指定序列类型：

```powershell
# 氨基酸序列：后续执行密码子优化
uv run python main.py add-auxiliary-protein -i config.json `
  --protein-file helper.txt --sequence-type protein

# CDS：后续直接使用，不执行密码子优化
uv run python main.py add-auxiliary-protein -i config.json `
  --protein-file helper_cds.fasta --sequence-type cds
```

FASTA 可以包含多条记录，ID 取 header 的第一个字段；纯文本按一条序列处理，ID 取文件名。多次导入会累积，同 ID 以最后一次上传的类型和内容为准。系统只去除空白并转为大写，不检查上传 CDS 的字符、三联体、起止密码子或内部终止密码子。

## 输入规则

- 主酶来自 `main_enzyme_selection.proteins`。
- 旧研究流程仍可从 `main_enzymes[*].confirmed_auxiliary_proteins` 读取辅助蛋白。
- 手动上传流程从 `auxiliary_protein_selection.proteins` 读取路线级共享辅助蛋白；这些蛋白不绑定具体主酶。
- `auxiliary_protein_selection.can_advance` 必须为 `true`；该区段不存在时只处理主酶并写入警告。
- 手动氨基酸序列从项目快照读取，不访问 UniProt；手动 CDS 直接写入 `cds_selection` 并标记 `optimization_skipped: true`。
- 目前仅支持 `ecoli_mg1655`，映射到 CodonTransformer organism ID 52。
- 主酶写入 manifest 时不再记录氨基酸序列哈希；运行本流程时直接读取或下载 FASTA。

## 下载、缓存与输出

输出目录为：

```text
outputs/<target>/protein_to_cds/
├── uploaded_sequences/manifest_revision_*/<id>.<type>.fasta
├── protein_sequences/<accession>.fasta
├── raw_cds/<accession>.raw.fasta
├── optimized_cds/<accession>.optimized.fasta
├── reports/<accession>.optimization.json
└── run_summary.json
```

如果本地蛋白 FASTA 是合法、单条且 accession 一致的记录，流程会直接复用，不访问 UniProt。无效或 accession 冲突的文件不会被覆盖。通过约束门禁且输入指纹一致的优化结果也会直接复用。

由氨基酸生成的 CDS 必须保持蛋白翻译完全一致，并通过起止密码子、全局及局部 GC、CAI、稀有密码子簇、双链禁用 motif 和同聚物检查。用户直接上传的 CDS 不经过这些优化门禁；现有表达盒和组装阶段仍保留自身的序列读取与完整性检查。

## manifest 与退出状态

运行结果会整体写入 `cds_selection`：

- `complete`：全部成功，CLI 退出码为 0。
- `partial`：部分成功；成功产物和失败原因都会写入 manifest，CLI 退出码为 2。
- `failed`：全部失败；所有失败原因都会写入 manifest，CLI 退出码为 2。

处理期间若 manifest revision 发生变化，本次结果不会覆盖新版本。重新选择路线、主酶或导入/替换辅助序列时，旧 `cds_selection` 会被自动清除。

# 复制模块库

`replication_modules.csv` 收录 BASIC SEVA 的5个复制模块。字段为 `id,name,host_range,copy_number,sequence`。

`sequence` 仅保存完整复制模块，不包含 T1、LMP/LMS、BSEVA_L1、抗性模块或 mScarlet 筛选盒。DNA 使用连续的大写 A/C/G/T 字符串，CSV 使用 UTF-8 编码。

- `host_range=escherichia_coli` 表示首版按大肠杆菌使用场景登记。
- `host_range=broad_host_range` 用于 RSF1010；其复制依赖模块内保留的 Rep 蛋白系统。
- `copy_number` 是面向选择的简化类别。实际拷贝数会随宿主、培养条件和完整质粒结构变化。
- `basic_seva_psc101_pkd46_ts` 是温敏版本：作者说明其在37°C不稳定，通常在30°C维持。

## 来源与提取范围

序列来自 BASIC SEVA 仓库提交 `cad60c243aa2b12c6ed52ce36be454da3f357a6e` 的 [BASIC_SEVA_collection_v10.gb](https://github.com/LondonBiofoundry/basicsynbio/blob/cad60c243aa2b12c6ed52ce36be454da3f357a6e/basicsynbio/parts_linkers/BASIC_SEVA_collection_v10.gb)。模块结构和温敏条件见 [BASIC SEVA 论文](https://doi.org/10.1093/synbio/ysac023)。

以下范围按来源记录的正链、1-based 闭区间表示：

| id | 来源记录 | 提取范围 | 长度（bp） |
| --- | --- | --- | ---: |
| basic_seva_rsf1010 | BASIC_SEVA_15a.10 | 1308–4978 | 3671 |
| basic_seva_p15a | BASIC_SEVA_16.10 | 1305–2036 | 732 |
| basic_seva_psc101 | BASIC_SEVA_17.10 | 1305–2764 | 1460 |
| basic_seva_pbr322_rop | BASIC_SEVA_19.10 | 1305–2689 | 1385 |
| basic_seva_psc101_pkd46_ts | BASIC_SEVA_17_pKD46.10 | 1297–2847 | 1551 |

pKD46 在主集合中以 `Rep101(pKD46)`、`pSC101 ori` 等子特征注释；该行保存 BSEVA_L1 结束后至 T1 开始前的整个复制模块。其边界还使用作者公开的 [独立核心序列](https://github.com/LondonBiofoundry/basicsynbio/blob/cad60c243aa2b12c6ed52ce36be454da3f357a6e/sequences/genbank_files/BASIC_SEVA_collection/misc_seqs/BS_x7x_pKD46_core_seq.gb) 交叉核对。

## 文件校验

已检查字段、行数、唯一编号、DNA 字符、模块边界、T1 排除情况及全部序列与来源片段的一致性。

- CSV SHA-256：`3d731f57397367dd17af25e5c8a64336e500c7d81eef39bb0c7356dd5dd074da`
- 来源 GenBank SHA-256：`a4cdb1eb8917f65f3727ec54108b604c97caa9768f546b9cf2ed3d2d9fc29917`

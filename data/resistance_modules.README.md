# 抗性模块库

`resistance_modules.csv` 收录 BASIC SEVA 的 6 类抗性、7 个模块版本。字段仅包含 `id,name,antibiotic,resistance_gene,sequence`。

`sequence` 是从原始 GenBank 中提取的连续片段：作者标注的完整抗性功能模块。原七条记录共同的 103 bp T0 已拆入 `terminator_modules.csv`；另将功能模块前部有来源依据的26/34 bp边界片段归入gap_modules.csv；保留功能模块内部调控与终止序列，不包含外围 LMP/LMS 接头、BSEVA_L1、复制模块或 mScarlet 筛选盒。DNA 使用连续的大写 A/C/G/T 字符串，CSV 使用 UTF-8 编码。加载时不会再按长度裁剪。

`antibiotic` 使用英文小写名称；链霉素／壮观霉素共用一个模块，写为 `streptomycin;spectinomycin`。

## 来源与提取范围

序列来源仓库固定到提交 `cad60c243aa2b12c6ed52ce36be454da3f357a6e`：

- [BASIC_SEVA_collection_v10.gb](https://github.com/LondonBiofoundry/basicsynbio/blob/cad60c243aa2b12c6ed52ce36be454da3f357a6e/basicsynbio/parts_linkers/BASIC_SEVA_collection_v10.gb)：30 个原始设计记录。
- [BASIC_SEVA_66-11_&_69-11.gb](https://github.com/LondonBiofoundry/basicsynbio/blob/cad60c243aa2b12c6ed52ce36be454da3f357a6e/basicsynbio/parts_linkers/BASIC_SEVA_66-11_%26_69-11.gb)：作者公布的庆大霉素变体记录。
- [BASIC SEVA 论文](https://doi.org/10.1093/synbio/ysac023)：组件结构、构建验证和版本条件。
- [原始 SEVA 论文](https://doi.org/10.1093/nar/gks1119)：抗性基因名称依据；四环素序列采用 BASIC SEVA 的 5a 版本。

以下范围按来源记录的正链、1-based 闭区间表示：

| id | 来源记录 | 提取范围 | 长度（bp） |
| --- | --- | --- | ---: |
| basic_seva_ap | BASIC_SEVA_15a.10 | 193–1231 | 1039 |
| basic_seva_km | BASIC_SEVA_25a.10 | 193–1119 | 927 |
| basic_seva_cm | BASIC_SEVA_35a.10 | 193–975 | 783 |
| basic_seva_sm_sp | BASIC_SEVA_45a.10 | 201–1189 | 989 |
| basic_seva_tet_5a | BASIC_SEVA_5a5a.10 | 201–1467 | 1267 |
| basic_seva_gm | BASIC_SEVA_65a.10 | 193–997 | 805 |
| basic_seva_gm_11 | BASIC_SEVA_66.11 | 191–1014 | 824 |

## 庆大霉素版本

`basic_seva_gm` 保留常规序列；`basic_seva_gm_11` 保留作者在 p15A/pBR322 组合中公布的序列变体。66.11 和 69.11 的该片段完全相同，因此合并为一条记录。不能据此假定两个版本与任意复制模块的组合均已验证；后续组装入口需处理版本匹配。

## 文件校验

已检查 CSV 的字段、行数、唯一编号、DNA 字符、模块边界及全部序列与来源片段的一致性。

CSV SHA-256：`b91ea292510cb9392e4217e42d20fe47a0f22413aae0f761410e509bd618cca4`。

拆分前 CSV SHA-256：`c6aed8ee3a2992975160dcc8d75598449145961d7f4c1fef9e0c3ab4b73327fb`。

来源文件 SHA-256：

- `BASIC_SEVA_collection_v10.gb`：`a4cdb1eb8917f65f3727ec54108b604c97caa9768f546b9cf2ed3d2d9fc29917`
- `BASIC_SEVA_66-11_&_69-11.gb`：`6d5e5b1cafbc688a0f4fa3cf2a9eb16482bf39b89afd9970022ae11e4dade90f`

Gap拆分原始快照、删除范围与坐标映射见src/plasmid_design/source_features.json及docs/gap-extraction-audit.md。

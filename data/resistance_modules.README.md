# 抗性模块库

`resistance_modules.csv` 收录 BASIC SEVA 的 6 类抗性、7 个模块版本。字段仅包含 `id,name,antibiotic,resistance_gene,sequence`。

`sequence` 是从原始 GenBank 中提取的连续片段：T0 终止子起点至完整抗性模块末端。保留中间原有调控和间隔序列，不包含外围 LMP/LMS 接头、BSEVA_L1、复制模块或 mScarlet 筛选盒。DNA 使用连续的大写 A/C/G/T 字符串，CSV 使用 UTF-8 编码。

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
| basic_seva_ap | BASIC_SEVA_15a.10 | 64–1231 | 1168 |
| basic_seva_km | BASIC_SEVA_25a.10 | 64–1119 | 1056 |
| basic_seva_cm | BASIC_SEVA_35a.10 | 64–975 | 912 |
| basic_seva_sm_sp | BASIC_SEVA_45a.10 | 64–1189 | 1126 |
| basic_seva_tet_5a | BASIC_SEVA_5a5a.10 | 64–1467 | 1404 |
| basic_seva_gm | BASIC_SEVA_65a.10 | 64–997 | 934 |
| basic_seva_gm_11 | BASIC_SEVA_66.11 | 62–1014 | 953 |

## 庆大霉素版本

`basic_seva_gm` 保留常规序列；`basic_seva_gm_11` 保留作者在 p15A/pBR322 组合中公布的序列变体。66.11 和 69.11 的该片段完全相同，因此合并为一条记录。不能据此假定两个版本与任意复制模块的组合均已验证；后续组装入口需处理版本匹配。

## 文件校验

已检查 CSV 的字段、行数、唯一编号、DNA 字符、模块边界及全部序列与来源片段的一致性。

CSV SHA-256：`c6aed8ee3a2992975160dcc8d75598449145961d7f4c1fef9e0c3ab4b73327fb`。

来源文件 SHA-256：

- `BASIC_SEVA_collection_v10.gb`：`a4cdb1eb8917f65f3727ec54108b604c97caa9768f546b9cf2ed3d2d9fc29917`
- `BASIC_SEVA_66-11_&_69-11.gb`：`6d5e5b1cafbc688a0f4fa3cf2a9eb16482bf39b89afd9970022ae11e4dade90f`

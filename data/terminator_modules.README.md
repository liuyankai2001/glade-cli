# 终止子模块库

`terminator_modules.csv` 仅包含 `id,name,sequence`，保存连续的大写 A/C/G/T 序列。

| id | name | 长度 bp | 来源记录 | 来源范围（1-based 闭区间） |
| --- | --- | ---: | --- | --- |
| basic_seva_t0 | T0 | 103 | BASIC_SEVA_15a.10 | 64–166 |
| basic_seva_t1 | T1 | 105 | BASIC_SEVA_15a.10 | 4993–5097 |

来源固定为 [BASIC SEVA 集合](https://github.com/LondonBiofoundry/basicsynbio/blob/cad60c243aa2b12c6ed52ce36be454da3f357a6e/basicsynbio/parts_linkers/BASIC_SEVA_collection_v10.gb)，提交 `cad60c243aa2b12c6ed52ce36be454da3f357a6e`。源文件 SHA-256 为 `a4cdb1eb8917f65f3727ec54108b604c97caa9768f546b9cf2ed3d2d9fc29917`。

T0 来自旧抗性模块共同的前 103 bp；T1 来自旧骨架元数据。角色、来源方向、功能注释和校验信息保存在 `src/plasmid_design/source_features.json`；CSV 是运行时终止子 DNA 的唯一来源。

用户需在网页中手动加入。参考排列为 `T1 → 完整表达构建 → T0`；缺项或错位会弹出警告，可继续生成当前序列。现有酶切冲突等阻断检查保持有效。

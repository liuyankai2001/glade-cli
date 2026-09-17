# Gap 来源及提取审计

范围：当前 7 个抗性、5 个复制、2 个终止子和固定接口；不扫描完整表达构建。

只拆出有来源依据且不覆盖功能注释的边界片段。未确认的非编码/调控区保留在模块中。所有原始 DNA、哈希、注释、删除范围、保留坐标映射及待核查区域存于 source_features.json。

| 模块 | 原始 bp | 当前 bp | 已拆出 bp | 保留待核查区域数 |
| --- | ---: | ---: | --- | ---: |
| basic_seva_ap | 1065 | 1039 | 26 | 1 |
| basic_seva_km | 953 | 927 | 26 | 1 |
| basic_seva_cm | 809 | 783 | 26 | 1 |
| basic_seva_sm_sp | 1023 | 989 | 34 | 1 |
| basic_seva_tet_5a | 1301 | 1267 | 34 | 1 |
| basic_seva_gm | 831 | 805 | 26 | 1 |
| basic_seva_gm_11 | 850 | 824 | 26 | 1 |
| basic_seva_rsf1010 | 3671 | 3671 | 0 | 0 |
| basic_seva_p15a | 732 | 732 | 0 | 0 |
| basic_seva_psc101 | 1460 | 1460 | 0 | 0 |
| basic_seva_pbr322_rop | 1385 | 1385 | 0 | 0 |
| basic_seva_psc101_pkd46_ts | 1551 | 1537 | 14 | 3 |
| basic_seva_t0 | 103 | 103 | 0 | 0 |
| basic_seva_t1 | 105 | 105 | 0 | 0 |

抗性前部26/34 bp依据作者的完整功能模块注释边界；pKD46尾部14 bp在集合及独立核心来源中紧邻SEVA_T1，且与原生复制接口完全一致。Ap末端含终止作用的DNA仍属于功能模块，未拆除。其他复制模块及T0/T1没有可确认的内部gap。pKD46局部1–61、1013–1059、1283–1537 bp可能参与调控，保留待核查。

归库 12 条完整DNA去重序列：BASIC L1–L6、BSEVA_L1、原生76/14 bp接口、24 bp占位及抗性26/34 bp接口。53 bp linker不含制备用GG接头；占位片段不自动加入最终质粒。相同DNA的全部本次提取来源保存在sources数组。

来源：[作者固定版本仓库](https://github.com/LondonBiofoundry/basicsynbio/tree/cad60c243aa2b12c6ed52ce36be454da3f357a6e)；[BASIC SEVA论文](https://doi.org/10.1093/synbio/ysac023)。原始文件保存在data/gap_sources，逐文件SHA256校验；重现：`python -B -m data.curate_gap_modules`。来源定义的中性linker是候选条目，不表示新位置组合已经实验验证。

# KEGG 全量筛选：尚未检索到大肠杆菌生成报道的候选目标

## 结论摘要

- 检索与计算日期：2026-08-21。
- KEGG 快照包含 19,623 个 Compound、12,473 个 Reaction 和 52,694 条 compound-reaction 原始链接。
- iML1515 底盘的 A0 包含 595 个唯一 KEGG 化合物。按照“所有反应底物必须已经可用”的规则扩展到 10 层后，累计可达 3,829 个化合物。
- 全量自动分档为：595 个 A0、3,234 个 1–10 层新增可达、5,741 个有 KEGG 反应但 10 层内不可达，以及 10,053 个没有直接 KEGG Reaction 的化合物。
- 严格人工核查后，保留 12 个“目标物尚未检索到大肠杆菌生成报道、且当前 GLADE 能返回路线”的 A 档候选。
- 另保留 13 个宿主新颖性较强、但当前 GLADE 被 KEGG 数据或搜索复杂度阻塞的 B 档候选。
- 没有候选同时达到“完整酶证据 + strict_l1 GEM 通过”的 A1 等级。11 个候选只能在 relaxed 辅因子模式产生通量；aphidicolin 虽然 strict_l1 有通量，但末端 KEGG 反应缺少酶锚点，不能视为可构建路线。

完整的 19,623 行筛选结果见 [ecoli_unreported_kegg_candidates.csv](ecoli_unreported_kegg_candidates.csv)。

## “未报道”的严格定义

以下任一情况都会判定目标物已经在大肠杆菌体系中生成，并从候选中排除：

1. 工程化大肠杆菌整细胞生成目标物；
2. 大肠杆菌无细胞或裂解液体系生成目标物；
3. 在大肠杆菌中表达的重组酶，经纯化或使用粗酶液后生成目标物。

仅使用大肠杆菌扩增质粒、表达只生成上游中间体的酶、以目标物作为底物生成衍生物，或用目标物进行抑菌实验，不算“大肠杆菌生成目标物”。

负面结论不能证明绝对不存在。因此本文统一使用：

> 截至 2026-08-21，在 PubMed、CrossRef/OpenAlex 及可访问全文中，尚未检索到大肠杆菌整细胞、无细胞体系或大肠杆菌表达酶生成目标产物的报告。

## A 档：当前 GLADE 可以搜索

### A2：GLADE 有路线，但仅 relaxed GEM 可行

| 优先级 | 目标物 | KEGG | 最短步数 | GLADE 解数 | strict_l1 | relaxed | 文献置信度 | 主要风险 |
|---:|---|---|---:|---:|---|---|---|---|
| 1 | Capsidiol | C09627 | 2 | 1 | FAIL | PASS | 高 | 植物 P450/CPR；通用电子载体未闭合 |
| 2 | Solavetivone | C09737 | 3 | 1 | FAIL | PASS | 高 | 同一 P450 连续两步氧化；膜表达与 CPR 风险 |
| 3 | Pisiferic acid | C09163 | 6 | 1 | FAIL | PASS | 高 | CYP76AH/CYP76AK 分支与产物特异性 |
| 4 | Salviol | C21819 | 6 | 1 | FAIL | PASS | 高 | P450 分支、副产物和 CPR 匹配 |
| 5 | Sugiol | C21822 | 6 | 2 | FAIL | PASS | 中 | CYP76AH 多产物反应 |
| 6 | 11-Hydroxyferruginol | C21796 | 6 | 2 | FAIL | PASS | 中 | 多功能 P450，容易继续氧化 |
| 7 | 11-Hydroxysugiol | C21823 | 7 | 2 | FAIL | PASS | 中 | 连续氧化和副产物风险 |
| 8 | 11,20-Dihydroxyferruginol | C21830 | 7 | 2 | FAIL | PASS | 中 | 两级 P450 氧化，电子传递未闭合 |
| 9 | Carnosic acid | C21818 | 7 | 5 | FAIL | 部分 PASS | 高 | 多个植物 P450；5 条路线中只有部分 relaxed 可行 |
| 10 | Oleanolic acid | C17148 | 7 | 1 | FAIL | PASS | 高 | 三萜膜酶、P450 与宿主表达负担高 |
| 11 | 11,20-Dihydroxysugiol | C21824 | 8 | 2 | FAIL | PASS | 中 | 最长的同族路线，多分支连续氧化 |

这里的 `PASS` 只代表在指定 GEM 假设下可以产生非零目标通量，不代表酶一定表达、底物特异性正确或实验一定成功。

### A3：能搜到，但反应证据不足

| 目标物 | KEGG | 最短步数 | strict_l1 | 问题 |
|---|---|---:|---|---|
| Aphidicolin | C06088 | 4 | PASS | 最后一步使用 R06316 或 R06317；KEGG 将这些反应标为 `unclear reaction`，并且缺少可用于蛋白选择的 EC/KO。当前 GLADE 输出却把路线标为 `resolved`，strict_l1 通过只是因为反应没有记录真实辅因子负担。 |

Aphidicolin 可以作为“发现数据库/解析器缺陷”的计算测试案例，不能直接进入蛋白选择和湿实验设计。

## 代表性路线

### Capsidiol：最短候选

```text
C00448 Farnesyl diphosphate
  --R09574 / EC 4.2.3.61-->
C19708 5-Epiaristolochene
  --R09573 / EC 1.14.14.149-->
C09627 Capsidiol
```

目标只需两个异源主酶，但第二步是植物 P450 的连续羟化反应，因此仍需解决 CPR、血红素和膜蛋白表达问题。已有工作在植物和酵母相关体系中验证该通路；本次严格检索未发现大肠杆菌体系形成 capsidiol 目标物的报告。

### Solavetivone：三步连续氧化候选

```text
C00448 Farnesyl diphosphate
  --R06523 / EC 4.2.3.21-->
C12142 Vetispiradiene
  --R09576 / EC 1.14.14.151-->
C19711 Solavetivol
  --R09577 / EC 1.14.14.151-->
C09737 Solavetivone
```

R09576 和 R09577 由同一 premnaspirodiene oxygenase 连续完成，适合测试 GLADE 能否把多个反应步骤归并到同一主酶，同时识别 CPR 需求。

### Phenolic abietane 路线家族

Pisiferic acid、salviol、sugiol、11-hydroxyferruginol、11-hydroxysugiol、两种 11,20-dihydroxy 产物和 carnosic acid 共享以下上游模块：

```text
GGPP -> Copalyl diphosphate -> Miltiradiene
     -> Abietatriene -> Ferruginol -> 不同 P450 氧化分支
```

这些不是 8 条完全独立的生物学路线，而是一个可复用上游模块对应多个末端氧化方案。它们适合做酶特异性和多方案表达设计比较，不应被当作 8 个互不相关的实验系统。

## B 档：宿主新颖性较强，但当前 GLADE 跑不通

| 目标物 | KEGG | 分类 | 当前阻塞原因 | 文献置信度 |
|---|---|---|---|---|
| Carnosol | C09069 | B1 | KEGG 没有目标 Reaction；涉及 carnosic acid 后续自发氧化 | 中 |
| Forskolin | C09076 | B1 | KEGG 没有目标 Reaction；大肠杆菌只做到 13R-manoyl oxide 等前体 | 高 |
| Betulinic acid | C08619 | B1 | KEGG 没有目标 Reaction；完整异源生产集中在酵母 | 高 |
| Ursolic acid | C08988 | B1 | KEGG 没有目标 Reaction | 中 |
| Parthenolide | C07609 | B1 | KEGG 没有目标 Reaction；大肠杆菌止于 costunolide，目标生产在酵母 | 高 |
| Withaferin A | C08841 | B1 | KEGG 没有目标 Reaction；完整下游酶网络仍不适合当前搜索 | 中 |
| Andrographolide | C20214 | B1 | KEGG 没有目标 Reaction；大肠杆菌只验证早期 CPP/copalol 模块 | 中 |
| Salvinorin A | C20196 | B1 | KEGG 没有目标 Reaction；通路仍在补全 | 中 |
| Artemisinin | C09538 | B1 | KEGG 没有目标 Reaction；大肠杆菌生产的是 amorphadiene/artemisinic acid 前体 | 高 |
| Quillaic acid | C08972 | B2 | R13150/R13151 标记为 `unclear reaction` | 高 |
| Ginkgolide A | C07601 | B2 | GLADE 实际拒绝 R08402：`unresolved_multistep_without_enzyme_anchor` | 高 |
| Colchicine | C07592 | B2 | 仅有末步 R08453，EC 为不完整的 `2.3.1.-`，上游仍缺失 | 高 |
| Paclitaxel | C07394 | B3 | KEGG 只有部分末步；批次进入该目标后直到 5 分钟总超时仍未完成，完整上游网络不闭合 | 高 |

其中 salvinorin A 具有精神活性，paclitaxel、withaferin A 和 colchicine 具有显著细胞毒性。这些条目仅作为计算与数据库覆盖案例，不构成实验实施建议。

## 严格定义下被排除的典型目标

| 目标物 | KEGG | 排除证据 |
|---|---|---|
| Momilactone A | C18015 | 在大肠杆菌表达的 momilactone-A synthase 酶测定中生成目标物 |
| Sclareol | C09183 | 已有工程化大肠杆菌整细胞生产报道 |
| Artemisinic acid | C20309 | 已有工程化大肠杆菌生产报道 |
| Glycyrrhetinic acid | C02283 | 已有大肠杆菌 whole-cell bioconversion 生成报道 |
| Linamarin | C01594 | 大肠杆菌表达的 UGT 裂解液酶反应生成目标物 |
| (+)-Camphor | C00808 | 大肠杆菌表达的 borneol dehydrogenase 酶反应生成目标物 |
| Noscapine | C09592 | 大肠杆菌表达并纯化的末步酶生成目标物 |
| Catharanthine | C09107 | 大肠杆菌表达的重组酶参与体外体系形成目标物 |

这也说明旧文档中按“没有全细胞从头生产”筛出的 camphor、noscapine 和 catharanthine，不满足本报告采用的更严格标准。

## 数据与方法边界

- 全量 CSV 对每个 KEGG Compound 给出自动状态：宿主 A0、1–10 层可达、有反应但 10 层内不可达、没有直接反应。
- `NOT_MANUALLY_AUDITED` 不表示未报道，只表示没有进入本轮人工全文核查短名单。
- A 档路线均由当前仓库代码实际运行生成，不是根据 KEGG 图手工猜测。
- GEM 验证统一使用 iML1515、默认培养基、per-solution 模式；为控制批量运行时间跳过 FVA，因此状态写作 `PASS_PRODUCT_FLUX_FVA_NOT_RUN`。
- relaxed 模式会开放通用电子载体，只能作为“补全辅因子后可能可行”的证据。实验设计前仍需显式加入 CPR、ferredoxin 或其他再生模块。
- PubMed 的部分批量请求出现 SSL EOF/握手超时，已通过重试、OpenAlex/网页搜索和出版社全文补查；因此负面检索仍可能遗漏未索引论文、学位论文、会议材料和未公开专利。

## 主要证据来源

- Carnosic acid 路线及酵母重建：[Nature Communications 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5059481/)
- Pisiferic acid 与 salviol 的定向酵母生产：[Scientific Reports 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5562805/)
- Tanshinone/abietane P450 分支：[New Phytologist 2016](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.13790)
- Solavetivone oxygenase：[Journal of Biological Chemistry 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC2695360/)
- Aphidicolin 在 Aspergillus oryzae 的完整异源重建：[Bioscience, Biotechnology, and Biochemistry 2011](https://www.jstage.jst.go.jp/article/bbb/advpub/0/advpub_110366/_pdf/-char/en)
- Forskolin 完整通路与大肠杆菌前体边界：[eLife 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5388535/)
- Parthenolide 的酵母生产及大肠杆菌 costunolide 边界：[Bioresources and Bioprocessing 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12141831/)
- Paclitaxel 的大肠杆菌早期通路：[Science 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC3034138/)
- Colchicine 前体在大肠杆菌中的人工级联：[Nature Communications 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10761944/)
- Momilactone A 的严格排除证据：[Journal of Biological Chemistry 2007](https://doi.org/10.1074/jbc.M703344200)

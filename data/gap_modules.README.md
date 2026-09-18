# Gap 序列库

CSV字段严格为id,name,sequence；DNA为大写ACGT。ID=gap_+完整序列SHA256；完全相同DNA合并，保留aliases及所有sources。

当前仅6条：BASIC L1–L6，均为53 bp通用中性linker候选，不含制备用GG接头。用途、证据及注释存于src/plasmid_design/source_features.json；模块原始快照及删除映射也存于该文件。已拆出的原生边界仅留存来源审计，不再作为可选Gap。详见docs/gap-extraction-audit.md。

新旧请求均只包含明确选入的gap实例，不自动补入旧76/14 bp间隔。临时骨架占位使用BASIC L1，最终由完整表达构建替换，不自动残留。纯序列拼接可重复选入同一条DNA；实际BASIC装配一轮内应使用不同linker。候选用于独立功能模块边界，不保证任意组合均具有生物学中性。重新运行归库脚本仍只生成L1–L6。

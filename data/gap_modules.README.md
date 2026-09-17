# Gap 序列库

CSV字段严格为id,name,sequence；DNA为大写ACGT。ID=gap_+完整序列SHA256；完全相同DNA合并，保留aliases及所有sources。

当前12条：6个BASIC标准linker、BSEVA_L1、76/14bp原生接口、24bp占位、26/34bp抗性接口。用途、证据及注释存于src/plasmid_design/source_features.json；模块原始快照及删除映射也存于该文件。详见docs/gap-extraction-audit.md。

新装配schema2只包含明确选入的gap实例；旧装配仅将原76/14bp外部间隔显式迁移，不恢复已经拆出的内部边界。可重复选入同一条DNA。增加数据时需同步固定来源、特征快照及用途；不要把任意非编码DNA推断为gap。

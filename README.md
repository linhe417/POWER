# Lithium solid-electrolyte conductivity dataset

公开实验数据的可追溯整合与清洗版本，2026-11-05。**五库全量尚未齐全；本项目尚未证实新物理规律。**

| 来源 | 主表来源记录 | 状态 |
|---|---:|---|
| OBELiX | 599 | 完整公开原表 |
| Liverpool LiIon | 820 | 完整官网表，含2024年修订 |
| McHaffie / Caltech | 571 | 完整公开原表 |
| Laskowski | 1,345 | 未核实二手转录；论文报告1,346条 |
| DDSE / DigBat | 0 | 官网限制公众整库复制，未获取 |

主表合计 **3,335条来源记录**，清洗子集 **1,878条精确数值记录**，室温候选 **1,393条**，实际CIF **321个**。来源之间存在重叠，记录数不等于独立样品数；同组成也不等于同一晶相或同一样品。

## 数据下载与浏览

- [完整来源主表 CSV](data/processed/unified/measurement_records.csv)
- [清洗后精确数值 CSV](data/processed/unified/clean_numeric_records.csv)
- [室温候选 CSV](data/processed/unified/reported_rt_candidates.csv)
- [Excel 工作簿](data/processed/unified/unified_electrolytes.xlsx)
- [321个CIF及其清单](data/processed/unified/structures/)
- [字段字典](data/processed/unified/data_dictionary.csv)
- [来源、许可和SHA256](data/processed/unified/source_manifest.json)
- [逐条清洗日志](data/processed/unified/cleaning_log.csv)
- [跨库候选对应](data/processed/unified/cross_source_links.csv)
- [完整数据说明](data/processed/unified/数据集说明.md)

主标签统一为S/cm；37条OBELiX删失界限与52条Laskowski转录界限均保留，未当作精确回归值。室温候选含来源仅声明RT的记录，温度、bulk/total、直接测量/外推等条件仍可能未知。公开CIF含坐标扰动。Laskowski未核实转录不会进入默认精确清洗子集。

## 离线复现

```bash
python3 -m pip install numpy pandas lxml openpyxl
python3 run_unified_dataset.py
python3 -m unittest discover -s tests -p test_unified_dataset.py -v
```

所需原始文件随仓库保存；复现不需要DDSE账号、ICSD访问、付费API、GPU或新的DFT计算。后续模型应使用`leakage_group_id`分组，避免共享组成或文献造成信息泄漏。

## 来源与使用条款

请阅读[DATA_LICENSES.md](DATA_LICENSES.md)。各数据源保留自身许可；本仓库不为整个混合集统一授予CC0或商业使用许可。

- OBELiX：[原仓库](https://github.com/NRC-Mila/OBELiX)，DOI [10.1039/D5DD00441A](https://doi.org/10.1039/D5DD00441A)，CC-BY-4.0。
- Liverpool：[官方完整表](https://pcwww.liverpool.ac.uk/~msd30/lmds/LiIonDatabase.html)，DOI [10.1038/s41524-022-00951-z](https://doi.org/10.1038/s41524-022-00951-z)。按官网要求引用University of Liverpool；商业使用需向该大学取得许可。
- McHaffie：[CaltechDATA](https://data.caltech.edu/records/23mvv-6gk43)，DOI [10.1039/D5DD00052A](https://doi.org/10.1039/D5DD00052A)，数据集CC0-1.0。
- Laskowski：DOI [10.1039/D2EE03499A](https://doi.org/10.1039/D2EE03499A)。本版采用OBELiX仓库附带转录，标记为未核实；未取得完整原始SI。
- DDSE：[官方平台](https://www.digbat.org/)。本仓库仅保留其获取状态，未收录整库数据。

此版本用于可追溯的数据整理和探索分析，尚未完成逐条原文复核，不构成已确认的新科学发现。

# Lasso 肽复合物结构证据补充（RCSB 直接核验）

> 通过 RCSB 检索 API 直接核验（2026-01），作为子代理调研的交叉验证。
> 检索词 "lasso peptide" 命中 80 个 PDB 条目，其中**与蛋白/核酸形成复合物**的条目如下。

## A. Lasso 肽 × 人类蛋白（直接正样本，最高价值）

| PDB | 复合物 | 方法/分辨率 | 文献 | 意义 |
|---|---|---|---|---|
| **9KDF** | 人内皮素 B 受体 ETB（calcineurin 融合）× **RES-701-3** | cryo-EM | Shihoya et al., Nat Commun 2025, [10.1038/s41467-025-57960-x](https://www.nature.com/articles/s41467-025-57960-x) | **首个 lasso 肽-人 GPCR 复合物**；RES-701-3 是 RES-701-1 同源物（ETB 选择性拮抗剂） |
| **9KDG** | 人 ETB 无配体形式（apo 对照） | cryo-EM | 同上 | 用于提取配体诱导接触差异/界面残基 |
| **7VF3** | 人 **Plexin B1** 胞外域 × lasso 骨架嫁接肽 **PB1m7** | X 射线 | *De novo Fc-based receptor dimerizers differentially modulate PlexinB1 function*, Structure 2022, [10.1016/j.str.2022.07.008](https://doi.org/10.1016/j.str.2022.07.008) | **lasso 拓扑肽-人受体**第二例（工程化 lasso 嫁接肽，PLXNB1 是抗癌/神经再生靶点） |
| **7VG7** | 人 Plexin B1 × lasso 嫁接肽 **PB1m6A9** | X 射线 | 同上 | 同靶点第二变体，可做变异体-亲和力配对样本 |

## B. Lasso 肽 × 细菌靶标（同源迁移依据）

| PDB | 复合物 | 方法/分辨率 | 文献 | 人同源迁移 |
|---|---|---|---|---|
| **6N60** | E. coli RNAP σ70 全酶 × **MccJ25** | cryo-EM | Braffman et al., PNAS 2019, PMID 30626643 | → POLR2A/POLR2B（RNAP II，次级通道+桥螺旋保守；α-amanitin 复合物 1K83/3CQZ/6EXV 作锚点） |
| **6N61** | E. coli RNAP σ70 全酶 × **Capistruin** | cryo-EM | 同上 | 同上 |
| **9DFC / 9DFD / 9DFE** | T. thermophilus **70S 核糖体** × **lariocidin**（3 个状态） | X 射线 2.5–2.6 Å | "A Broad Spectrum Lasso Peptide Antibiotic Targeting the Bacterial Ribosome", Res Sq 2024, [10.21203/rs.3.rs-5058118/v1](https://doi.org/10.21203/rs.3.rs-5058118/v1) | → 人 80S 核糖体（肽基转移酶中心/A 位 23S↔28S rRNA 高度保守）；**新迁移对** |
| **4CU4** | E. coli 外膜转运蛋白 **FhuA** × MccJ25 | X 射线 | MccJ25 摄取途径（铁载体转运体） | → 人无直接同源（细菌外膜特异），可作"肽-转运蛋白"弱样本 |

## C. 生物合成/成熟酶 × 前体肽（可用于"酶-底物"迁移、注意力对齐弱监督）

| PDB | 复合物 | 说明 |
|---|---|---|
| 5V1U / 5V1V | TbiB1 × TbiA(β/α) leader peptide | 异肽键合成酶-前体复合物 |
| 6JX3 | Lasso 肽合成酶 B1 × leader peptide | 同源复合物 |
| 9X90 | PbaB1 × PbaA leader peptide | 同源复合物 |
| 8ITG / 8GQA | 差向异构酶 MslH × 前体肽 MslA 变体 | 修饰酶-底物 |
| 5TXE | 异肽酶 AtxE2(S527A) × Astexin3 | 水解酶-产物复合物（逆反应） |
| 8IBO | Mtb ClpC1 NTD × Lassomycin（子代理报告已含） | 靶标复合物 |

## D. 与 LassoPred 数据库的衔接

- LassoPred 数据库（4749 条）中 28 条带 PDB_ID；上述复合物涉及的肽序列可直接回表 `lassopred_database.csv` 的 `PDB_ID` 列。
- **9KDF 中的 RES-701-3 即库内条目 LP_434**（core: `GNWHGTSPDWFFNYYW`）——"数据库条目 ↔ 复合物结构"直接打通。
- 本地库内其他已知活性肽检索结果：anantin 4 条（LP_366/371/387 等）、siamycin I/III（LP_37/LP_43，后者 PDB 1RPB）、capistruin 3 条、microcin J25（LP_471，PDB 1Q71）、chaxapeptin（LP_333，PDB 2N5C）、sphingopyxin 11 条——构成"文献活性肽 × 库内序列"映射，便于组扫描集。

## E. 直接人类靶标结构正样本清单（汇总两条子代理调研后）

| 肽 | 人类靶标 | PDB | 证据强度 |
|---|---|---|---|
| RES-701-3 | ETB (EDNRB) | 9KDF | ★★★ 直接复合物 |
| PB1m7 / PB1m6A9（lasso 嫁接） | Plexin B1 (PLXNB1) | 7VF3 / 7VG7 | ★★★ 直接复合物 |
| α-amanitin（双环肽，非 lasso） | POLR2A（Pol II） | 6EXV / 3CQZ / 1K83 | ★★☆ 同源界面迁移锚点 |
| MccJ25 / Capistruin | 细菌 RpoC → 人 POLR2A | 6N60 / 6N61 | ★★☆ 同源迁移 |
| lariocidin | 细菌 23S rRNA → 人 28S rRNA | 9DFC/D/E | ★★☆ 同源迁移 |
| CsA / FK506 / compstatin 等（环肽，非 lasso） | PPIA / FKBP1A / C3 | 1CWA 等 | ★★☆ 跨拓扑迁移 |

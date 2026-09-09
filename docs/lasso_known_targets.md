# Lasso 肽已知靶标 / 生物活性调研汇总

> 用途：为「Lasso 肽–人类蛋白互作预测」项目构建正样本集。
> 调研方式：web_search 文献检索（2026 会话）；未使用本地全文库，序列/结合常数以公开文献与数据库条目为准，个别字段标注待核对。
> 证据等级：**A** = 有结构复合物（最高置信）；**B** = 生化/结合实验直接证据；**C** = 表型/活性（分子靶标未直接验证）。

---

## 1. 直接作用于人类蛋白的 Lasso 肽（最高优先级）

| 肽名 | 序列（如找到） | 靶标蛋白（物种 + UniProt） | 证据类型 | 复合物 PDB | 关键文献（PMID / DOI） |
|---|---|---|---|---|---|
| **RES-701-1** | `GNWHGTAPDWFFNYYW`（16 aa；环 Gly1–Asp9，C 端尾穿环） | 人内皮素 B 受体 ETB / EDNRB（**P24530**），GPCR；选择性拮抗剂（>1000× 相对 ETA） | 放射配体结合（B）+ 同源肽 RES-701-3 的冷冻电镜复合物（A） | **9KDF**（ETB–RES-701-3，EMDB EMD-62272/62273；注意复合物用的是同族肽 RES-701-3，非 RES-701-1 本身） | PMID 8175478（I 生产菌/发酵/性质）、8175479（II 序列）、8183252（选择性 ETB 拮抗剂）、8564420（溶液 NMR 结构）；Shihoya et al. 2025 *Nat Commun* DOI 10.1038/s41467-025-57960-x（ETB 复合物、GPCR 反向激动机制） |
| **Anantin** | 17 aa lasso 肽（一级序列见 PMID 1826288 / DRAMP18339 / IUPHAR ligand 4866；序列与环化方式待核对） | 人心钠素受体 A **NPR-A / NPR1 / GC-A**（**P16066**，鸟苷酸环化酶偶联受体）；ANP/ANF 竞争性拮抗剂 | 竞争性结合/功能拮抗（B） | 无 | PMID 1826287（I 生产菌/分离/活性）、1826288（II 序列，*J Antibiot* 1991 44:172-180，DOI 10.7164/antibiotics.44.172） |
| **Propeptin**（Propeptin-1） | 环肽（Microbispora 来源） | 人脯氨酰寡肽酶 / 脯氨酰内肽酶 **PREP**（**P48147**）；酶抑制剂 | 酶活抑制（B） | 无 | Esumi et al. 2002/2003（结构测定）；Kimura et al. 2007（Propeptin-2 类比物活性） |
| **Chaxapeptin** | lasso 肽（Streptomyces leeuwenhoekii C58） | 表型活性：抑制人肺癌 A549 细胞侵袭；**无明确分子靶标** | 表型（C） | 无 | PMID 26402731（*J Org Chem* 2015, 80:10252, DOI 10.1021/acs.joc.5b01878） |
| **Sungsanpin** | lasso 肽（深海 Streptomyces sp.） | 表型活性：抑制人肺癌 A549 细胞侵袭；**无明确分子靶标** | 表型（C） | 无 | PMID 23662937（*J Nat Prod* 2013） |
| **Ulleungdin** | lasso 肽（基因组挖掘，Streptomyces） | 表型活性：抑制癌细胞迁移；**无明确分子靶标** | 表型（C） | 无 | DOI 10.1021/acs.jnatprod.8b00449（*J Nat Prod* 2018） |
| **Sphingopyxin I** | lasso 肽（Sphingopyxis alaskensis） | 表型活性：抗癌/细胞毒（人癌细胞系）；**无明确分子靶标**（其异肽酶结构已解析，见 Nebel et al. *Angew Chem* 2016, DOI 10.1002/anie.201605232） | 表型（C） | 无 | Nebel et al. 2016 *Angew Chem Int Ed*（Sphingopyxin I 异肽酶结构） |

> 说明：RES-701-1 是当前唯一有 **人类蛋白复合物结构** 的 lasso 肽（经同族肽 RES-701-3 落地的 ETB 结构），是本项目最硬的直接人类靶标正样本。Anantin（NPR-A）与 Propeptin（PREP）为有直接结合/酶活证据的人类靶标，但缺复合物结构。抗癌类（Chaxapeptin/Sungsanpin/Ulleungdin/Sphingopyxin I）仅有表型活性，需在下游单独作为「线索级」样本处理。

---

## 2. 抗菌 Lasso 肽的细菌靶标（同源迁移依据）

| 肽名 | 序列（如找到） | 靶标蛋白（物种 + UniProt） | 证据类型 | 复合物 PDB | 关键文献（PMID / DOI） |
|---|---|---|---|---|---|
| **Microcin J25（MccJ25）** | `GAGHVPEYFVGIGTPISFYG`（21 aa；环 Gly1–Glu8） | 细菌 RNA 聚合酶 β′ 亚基 **RpoC** 次级通道（*E. coli* RpoC = **P0A8T7**）；阻塞次级通道抑制转录延伸、阻断 NTP 进入 | 生化（B）+ 冷冻电镜复合物（A） | **6N62**（*E. coli* RNAP–MccJ25，Braffman 2019）；肽单体：1PP5（X 射线）、1Q71（NMR） | PMID 15200952（Mukhopadhyay et al. *Mol Cell* 2004, DOI 10.1016/j.molcel.2004.06.010）、30626643（Braffman et al. *PNAS* 2019, DOI 10.1073/pnas.1817352116） |
| **Capistruin** | `GTPGFQTPDARVISRFGFN`（19 aa；环 Gly1–Asp9） | 细菌 RNA 聚合酶次级通道（*Burkholderia thailandensis* 产） | 生化（B）+ 冷冻电镜复合物（A） | **6N61**（*E. coli* RNAP–Capistruin）；肽单体 NMR：2KRN（BMRB 20014） | PMID 18616307（Knappe et al. *JACS* 2008，发现）；30626643（Braffman et al. *PNAS* 2019，复合物结构） |
| **Acinetodin** | lasso 肽 | 细菌 RNA 聚合酶（次级通道靶向） | 生化（B） | 无 | PMID 28106375（Metelev et al. *ACS Chem Biol* 2017, 12:814, DOI 10.1021/acschembio.6b01154） |
| **Klebsidin** | lasso 肽 | 细菌 RNA 聚合酶（次级通道靶向） | 生化（B） | 无 | PMID 28106375（同上；产自临床分离肺炎克雷伯菌 *Klebsiella pneumoniae*） |
| **Ubonodin** | lasso 肽 | 细菌 RNA 聚合酶；对洋葱伯克霍尔德复合群（*Burkholderia cepacia* complex）具抗菌活性 | 生化（B） | 无 | 发现：Cheung-Lee et al. *ACS Chem Biol* 2019；转运通路：PMID 35802499（*ACS Chem Biol* 2022, 17:2332, DOI 10.1021/acschembio.2c00420） |
| **Xanthomonin I/II** | lasso 肽（Xanthomonas 产） | 细菌 RNA 聚合酶（次级通道靶向） | 生化（B） | 无 | Hegemann et al.（Xanthomonas 基因组挖掘）；具体 PMID 待核对 |
| **Citrulassin** | lasso 肽（含瓜氨酸） | 报道为 RNA 聚合酶靶向（**待核对**） | 生化（B，待核对） | 无 | 待核对 |
| **Lassomycin** | lasso 肽（*Lentzea kentuckyensis*） | 结核分枝杆菌 **ClpC1 ATPase**（**P9WPC9**，N 端结构域）；刺激 ATP 酶、使 ATP 水解与蛋白降解解偶联 → 杀菌 | 生化（B）+ 抗结核表型（C） | 肽单体 NMR：**2MAI** | PMID 24684906（Gavrish et al. *Chem Biol* 2014, 21:509, DOI 10.1016/j.chembiol.2014.01.014） |
| **Lariatin A/B** | lasso 肽（*Rhodococcus* sp. K01-B0171） | 抗结核（*M. tuberculosis* 生长抑制）；**分子靶标未明确** | 表型（C） | 无 | Iwatsuki et al. 2006 *JACS*（Lariatin 具 lasso 结构）；序列见 peptidedb「Lariatin A」 |
| **Humidimycin（MDN-0010）** | lasso 肽（*Streptomyces humidus*） | 抗真菌：卡泊芬净（caspofungin）活性增强剂/协同剂（推测影响真菌细胞壁/葡聚糖合成通路）；**分子靶标未精确定义** | 表型（C） | 无 | Valiante et al. 2017（生物合成基因簇，PMC7168211） |

> 说明：RNA 聚合酶（RpoC 次级通道）是 lasso 肽最明确的细菌靶标族，且 MccJ25 与 Capistruin 已有复合物结构（6N62/6N61），可作为「肽–蛋白界面」的结构同源模板；Lassomycin 是唯一已知靶向 ClpC1 ATPase 的 lasso 肽。Lariatin/Humidimycin 等仅有表型活性。

---

## 3. 抗病毒 Lasso 肽（Siamycin 家族）

| 肽名 | 序列（如找到） | 靶标蛋白 | 证据类型 | 复合物 PDB | 关键文献（PMID / DOI） |
|---|---|---|---|---|---|
| **Siamycin I**（= RP 71955 = Aborycin = NP-06） | `CLGIGSCNDFAGCGYAIVCFW`（21 aa，2 对二硫键，三环） | HIV-1 包膜糖蛋白 **gp41 / gp120**（Env，HIV-1 HXB2 gp160 = **P04578**）；抑制合胞体形成/膜融合（病毒进入） | 表型（C）+ 溶液结构（A，肽自身） | 肽溶液结构：**1RPB** | Fréchet et al. 1994（RP 71955 溶液结构）；Tsunakawa et al. 1995（Siamycin I / NP-06 抗 HIV）；Detlefsen et al. 1995 |
| **Siamycin II** | `CLGVGSCNDFAGCGYAIVCFW`（21 aa，2 对二硫键；与 I 的第 4 位 I/V 差异） | HIV-1 包膜糖蛋白 **gp41/gp120**（融合抑制） | 表型（C）+ 溶液结构（A，肽自身） | 无复合物 | PMID 7787424（Constantine et al. *J Biomol NMR* 1995, DOI 10.1007/BF00211754，Siamycin II 溶液结构） |

> 说明：Siamycin I/II 的靶标为病毒（HIV）蛋白而非人类蛋白，但作用界面是人类细胞–病毒融合过程；现有证据为「抑制 HIV-1 合胞体形成 + 与 gp41/gp120 相互作用」的表型/生化线索，未见肽–gp41 复合物结构。

---

## 4. 其他实验验证活性（汇总）

| 肽名 | 活性 | 靶标/机制 | 证据等级 | 文献线索 |
|---|---|---|---|---|
| Propeptin | 酶抑制 | 人脯氨酰寡肽酶 PREP（见第 1 节） | B | Kimura/Esumi 系列 |
| Lassomycin | 抗结核 | ClpC1 ATPase（见第 2 节） | B | PMID 24684906 |
| Lariatin A/B | 抗分枝杆菌 | 未明确 | C | Iwatsuki 2006 *JACS* |
| Humidimycin | 抗真菌（卡泊芬净增效） | 未明确（细胞壁通路） | C | Valiante 2017 |
| Chaxapeptin / Sungsanpin / Ulleungdin / Sphingopyxin I | 抗癌（侵袭/迁移抑制、细胞毒） | 未明确 | C | 见第 1 节 |
| MS-271 | 抗革兰阳性菌（含 MRSA） | 未明确 | C | Yano et al. 1996（待核对） |
| Astexin-1/2/3、Caulonodin、Benenodin、Sphingonodin、Svicin 等 | 多无明确生物活性报道（多用于折叠/热稳定性/[1]轮烷开关研究） | — | — | Zong et al. 2017 *JACS*（Benenodin-1 热驱动 [1]轮烷开关） |

> 注：**抗疟活性**：本次未检索到有明确抗疟靶标/机制的 lasso 肽条目，列为发现缺口。

---

## 5. 综述推荐

| 综述 | 期刊/年份 | 说明 | 链接/ID |
|---|---|---|---|
| Maksimov, Pan & Link | *Nat Prod Rep* 2012 | Lasso 肽结构、功能、生物合成与工程（含经典活性表） | PMID 22833149 / DOI 10.1039/c2np20070h |
| Hegemann, Zimmermann, Xie, Marahiel | *Acc Chem Res* 2015 | 「Lasso Peptides: An Intriguing Class of Bacterial Natural Products」 | PMID 26079760 |
| Zong et al. | *JACS* 2017 | Lasso 肽折叠/[1]轮烷热开关（Benenodin-1） | DOI 10.1021/jacs.（待核对期卷页） |
| （Kodani 相关） | *J Antibiot* 等 | Streptomyces lasso 肽综述 | 待核对 |
| 「Put a Bow on It: Knotted Antibiotics Take Center Stage」 | *Antibiotics* 2019, 8(3):117 | Lasso 肽抗菌综述 | DOI 10.3390/antibiotics8030117（PMC6784204） |
| 「Lasso Peptides: Heterologous Production and Potential Medical Application」 | *Front Bioeng Biotechnol* 2020 | 含抗病毒等活性表（Table 3） | PMC7549694 / DOI 10.3389/fbioe.2020.571165 |

---

## 6. 可作为正样本的「直接人类靶标互作对」清单

以下为**肽 ↔ 人类蛋白**的直接互作对（按证据强度排序）：

1. **RES-701-1 ↔ EDNRB（P24530）** — 证据 A/B。复合物结构存在（经 RES-701-3 同族肽，PDB 9KDF）。**首选正样本**。
2. **RES-701-3 ↔ EDNRB（P24530）** — 证据 A。PDB 9KDF 直接沉积（calcineurin 融合人 ETB–RES-701-3）。与 RES-701-1 同族、可作增强正样本。
3. **Anantin ↔ NPR1 / NPR-A（P16066）** — 证据 B。ANP 竞争性拮抗。**正样本（无结构）**。
4. **Propeptin ↔ PREP（P48147）** — 证据 B。酶活抑制。**正样本（无结构）**。

以下为**表型线索级**（人类癌细胞活性，分子靶标未验证，暂不作为直接互作正样本，可作弱监督/迁移线索）：
- Chaxapeptin / Sungsanpin / Ulleungdin / Sphingopyxin I → 人肺癌/迁移/细胞毒表型。

以下为**同源迁移依据（非人类靶标）**，可用于验证「肽–蛋白界面」通用性：
- MccJ25 ↔ RpoC（P0A8T7，6N62）；Capistruin ↔ RNAP（6N61）；Lassomycin ↔ ClpC1（P9WPC9，2MAI 为肽单体）。

---

## 7. 发现缺口（后续补查建议）

1. **Anantin 精确序列与环化方式**：需核对 PMID 1826288 原文或 DRAMP18339 / IUPHAR 4866 的一级序列与 N 端 Gly–侧链羧基成环位点。
2. **RES-701-1 精确 Kd/IC50**：文献报道为「低 nM、对 ETB 高度选择性（相对 ETA >1000×）」，需补精确数值（结合常数来源：PMID 8183252 / 8175478）。
3. **抗癌 lasso 肽的分子靶标**：Chaxapeptin/Sungsanpin/Ulleungdin/Sphingopyxin I 目前仅有侵袭/迁移/细胞毒表型，无靶标——是「Lasso 肽–人类蛋白互作」预测最有价值的待填补区。
4. **Siamycin I/II ↔ gp41/gp120 的结合/结构证据**：仅有溶液结构 + 融合抑制表型，无肽–包膜糖蛋白复合物结构；是否属于「直接互作」需谨慎。
5. **Citrulassin、Xanthomonin 的靶标确认与序列**：本次证据不完整，需补原始文献。
6. **Lariatin/Humidimycin 的分子机制**：仅有抗结核/抗真菌表型，靶标未定。
7. **抗疟活性 lasso 肽**：本次未检索到明确条目。
8. **Zong 2017 JACS 精确卷期页**：待核对（本表仅标注 DOI 前缀）。

---

*文件生成于 2026 会话；所有未确认字段均已标注「待核对」，引用请以原始文献为准。*

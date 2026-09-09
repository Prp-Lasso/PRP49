# Lasso 肽–人类蛋白互作预测：可迁移同源样本调研报告

> 目标：为「Lasso 肽 ↔ 人类蛋白」互作预测模型寻找**可迁移的训练正样本**。
> 思路分两类：
> 1. **同源迁移**（homology transfer）：Lasso 肽已知结合的细菌/病毒靶标 → 具有结构同源关系的人类蛋白，重点是**结合界面**是否保守。
> 2. **跨拓扑迁移**（cross-topology transfer）：非 lasso 的环肽/大环肽药物与人类蛋白的**直接复合物结构**，作为「环肽–人靶」界面模式的正样本。
>
> 调研时间：以 web 检索到的公开结构数据库（PDB）与文献（PubMed/DOI）为准。

---

## 一、细菌 RNA 聚合酶 → 人 RNA 聚合酶 II（高价值）

### 1.1 关键事实

| 项目 | 内容 |
|---|---|
| 细菌靶标 | 大肠杆菌 RNAP **β′ 亚基 RpoC**（UniProt **P0A8T7**, RPOC_ECOLI）；MccJ25 主要与 β′ 的**次级通道（secondary channel）** 相互作用，capistruin 亦然 |
| 人同源物 | RNA Pol II 最大亚基 **POLR2A / RPB1**（UniProt **P24928**） |
| 结合界面 | β′/RPB1 的**次级通道（secondary channel / pore）**：位于活性中心附近、**桥螺旋（bridge helix）与触发环（trigger loop）之间**的窄通道，是 NTP 进入/产物排出的路径。α-amanitin 结合在**漏斗（funnel）内、桥螺旋下方**，主要由最大亚基残基构成 |
| 直接环肽证据 | **α-amanitin（鹅膏蕈碱）** 是双环八肽（bicyclic octapeptide），直接与人 POLR2A 复合物结构已解析 —— 这是「环肽–人 RNAP II」**直接复合物**证据 |

### 1.2 关键 PDB

- **6N60**：E. coli RNAP σ70 全酶 + 上游叉形启动子 DNA + **Microcin J25 (MccJ25)**（lasso 肽–RNAP 复合物）
- **6N61**：E. coli RNAP σ70 全酶 + **Capistruin**（lasso 肽–RNAP 复合物）
- **1K83**：α-amanitin–酵母 RNA Pol II 复合物（经典结构，Bushnell 等 2002）
- **3CQZ**：10 亚基 RNA Pol II + α-amanitin（Brueckner & Cramer 2008）
- 另见：哺乳动物 RNA Pol II 延伸复合物 + α-amanitin 的 cryo-EM 结构（JBC 2020）

### 1.3 空间关系（MccJ25 vs amanitin）

两者都作用在**双大亚基催化核心**上，且都与**桥螺旋（bridge helix）** 相邻，但结合口袋**相邻而非完全相同**：

- MccJ25/capistruin 的 lasso 环**插入次级通道**，堵住 NTP 进出路径；
- α-amanitin 结合在**漏斗（funnel）中、桥螺旋下方**的疏水口袋。
- 文献明确把二者并列：thermus RNAP 抑制剂 streptolydigin 与 MccJ25，以及酵母 RNAP II 抑制剂 α-amanitin，都**在桥螺旋附近**相互作用（PMID 20949106 / PMC2922424）。

### 1.4 序列保守性

- 细菌 β′（RpoC）与人 RPB1（POLR2A）**全长序列一致性较低（约 20–25%）**，但催化核心的「蟹钳（crab claw）」双 ψ–β 桶折叠保守；
- **桥螺旋、触发环、开关区（switch regions）及不变 DFDGD 基序**在细菌/古菌/真核三类生命中**高度保守**（触发环是「three domains of life 高度保守」的活性中心元件）。

### 1.5 迁移评级：**高**

> 桥螺旋 + 次级通道是跨域保守的活性中心地标，且已有 amanitin–人 POLR2A 的**直接环肽复合物**作为界面锚点；MccJ25/capistruin 与 amanitin 位点相邻、共享桥螺旋界面，迁移价值最高。

---

## 二、ClpC1 → 人 AAA+ ATPase（中价值）

### 2.1 关键事实

| 项目 | 内容 |
|---|---|
| 细菌靶标 | 结核分枝杆菌 **ClpC1**（UniProt **P9WPC3**，CLPC1_MYCTU），ClpC1P1P2 蛋白酶机器解折叠亚基 |
| 人同源物 | **CLPB**（caseinolytic peptidase B / 线粒体 disaggregase，UniProt **Q9H078**，别名 Skd3）；**CLPX**（UniProt **O76031**）；线粒体 ClpXP 机器 |
| 结合界面 | **Lassomycin 结合 ClpC1 的 N 端结构域（NTD）** —— 结构已解析（见下）。注意：人类的 CLPB NTD 是**锚蛋白重复（ankyrin repeat）**结构，CLPX NTD 是**锌结合结构域（ZBD）**，与 ClpC1 的螺旋状 NTD **拓扑不同**；保守的是 C 端 **AAA+ ATPase 核心（D1/D2 两个 AAA+ 域）** |
| 其他 AAA+ 靶向肽 | bortezomib 类（非肽天然产物也可作对照）、Bz-Leu-Leu 激活剂等（ClpC1P1P2 复合物结构见下） |

### 2.2 关键 PDB

- **8IBO**：野生型 M. tuberculosis **ClpC1 NTD + Lassomycin** 复合物（X 射线晶体结构）
- **8YD1**：M. tuberculosis ClpC1P1P2 复合物 + bortezomib（cryo-EM）
- **9IF4**：M. tuberculosis ClpC1P1P2 + 激活剂 Bz-Leu-Leu
- **7XC5**：人 **CLPB 的 ANK（锚蛋白）结构域**晶体结构
- **7US2 / 7TTR**：人 Skd3（CLPB）核苷酸结合域六聚体 / AAA+ 域

### 2.3 保守性评估

- AAA+ 解折叠酶核心（D1/D2）在 ClpC1 与人 CLPB/CLPX 之间**中等保守（约 30–40%）**；
- 但 **Lassomycin 的结合位点（NTD）在人类同源物中不保守**（螺旋 NTD vs 锚蛋白 / 锌结合域），因此「肽–NTD」界面**不能直接迁移**，只能迁移「肽–AAA+ 核心」层面的粗略模式。

### 2.4 迁移评级：**中**

> AAA+ ATPase 核心有同源性、可作为结构特征来源，但 Lassomycin 的 NTD 结合界面在人体无对应结构，界面级迁移受限。

---

## 三、脯氨酰寡肽酶（Propeptin 靶）→ 人 PREP（低价值，简略）

| 项目 | 内容 |
|---|---|
| 细菌靶标 | 细菌/放线菌（如 Microbispora）脯氨酰内肽酶（prolyl oligopeptidase, POP，S9 丝氨酸蛋白酶家族） |
| 人同源物 | **PREP**（脯氨酰内肽酶，UniProt **P48147**） |
| 结合界面 | POP 由 **β–螺旋桨 + α/β 水解酶催化域** 组成；Propeptin 为其抑制剂（文献以 Flavobacterium POP 作经典 assay 酶） |
| 关键 PDB | 人 PREP 结构如 **1QFM / 1H2W**（S9 家族）；细菌 POP 结构若干 |
| 文献 | Propeptin 发现（J. Antibiot., Kimura 等）；Prolyl endopeptidase inhibitors from actinomycetes（Biosci. Biotechnol. Biochem. 61:1754） |
| 评级 | **低**：S9 催化域跨物种保守，但 Propeptin–POP 的界面级复合物结构稀缺，迁移证据薄弱 |

---

## 四、跨拓扑迁移：环肽/大环肽药物–人类蛋白复合物（最高价值，重点）

> 这些是**非 lasso 环肽/大环肽直接结合人类蛋白**的高分辨率复合物，用于迁移学习的核心论据是：**拓扑（lasso 打结 vs 首尾环化 vs 双环）不同，但「肽主链–蛋白界面」的相互作用模式（疏水埋藏 + 氢键/盐桥 + 形状互补）是普适的**。

### 4.1 环孢素 A – 亲环蛋白 A（PPIA）

- **CsA**：11 残基环状十一肽（含 MeBmt 等非蛋白氨基酸）
- **PPIA / 亲环蛋白 A**：UniProt **P62937**
- **PDB**：**1CWA**（CsA–CypA）、**2RMA**（CypA + CsA 及衍生物，高分辨率）、**1CYH**
- **文献**：CsA–CypA 复合物晶体结构（Nature 1992 等）；CsA–CypA 复合物再结合钙调磷酸酶 calcineurin
- **价值**：环肽结合「环孢素结合口袋」的典型界面，正样本极丰富

### 4.2 FK506 / 雷帕霉素 – FKBP12

- **FKBP12 / FKBP1A**：UniProt **P62942**
- **PDB**：FK506–FKBP12 **1FKF / 1FKD / 1FKJ**；雷帕霉素–FKBP12 **1FAP / 1FKB**
- **价值**：大环内酯–脯氨酰异构酶界面的经典正样本；FKBP12 是**肽基脯氨酰顺反异构酶（PPIase）**，本身是肽结合蛋白

### 4.3 其他天然环肽 / 环肽–人靶复合物

| 肽 | 人靶蛋白 | PDB | 说明 |
|---|---|---|---|
| **Compstatin**（13 残基环肽） | 补体 **C3**（C3c） | **2QKI**（C3c–compstatin）、1A1P（NMR） | 环肽–人补体蛋白直接复合物，界面明确 |
| **MCoTI-II**（cyclotide, Möbius 族） | 胰蛋白酶（trypsin，牛源，人同源高度保守） | **4GUX** | 环肽 knottin 结合丝氨酸蛋白酶口袋 |
| **双环肽（bicyclic peptide）** | 转铁蛋白受体 **TfR1** | **9GH7** | 噬菌体展示双环肽–人 TfR1 |
| 双环肽 BCY00016132 | **NKp46**（NCR1） | **9H8R** | 双环肽–人受体复合物 |
| 双环肽抑制剂 | **ACE2** | **8BFW** | 双环肽–人 ACE2 |

### 4.4 讨论：能否作为迁移学习正样本？

- **支持（界面普适性）**：这些复合物证明「合成/天然环肽能以高亲和力结合人类蛋白的沟槽、口袋与环区」，其界面统计特征（接触残基类型、埋藏面积、极性互补）与 lasso 肽结合界面**同类可比**，可显著扩充「肽–人靶」正样本池。
- **保留（拓扑差异）**：lasso 肽的**打结（threaded）拓扑**赋予其「刚性 loop + 伸出的 tail」双模块结合方式，与首尾环化/双环肽的构象灵活性不同；因此建议在特征层显式编码**拓扑类别标签（lasso / cyclic / bicyclic / disulfide）**，让模型在共享界面特征的同时不混淆拓扑。
- **结论**：可作为**高置信度的跨拓扑正样本**，评级 **高**（界面迁移），但需配合拓扑标签做 domain adaptation。

---

## 五、病毒靶标迁移：Siamycin 抗 HIV（低价值）

| 项目 | 内容 |
|---|---|
| 病毒靶标 | HIV-1 **包膜糖蛋白 Env（gp160）**，融合机制（gp41/gp120 区）。UniProt 例：**P04578**（HIV-1 HXB2 Env） |
| 人同源物 | **无**（病毒 I 型跨膜融合蛋白，人类无同源蛋白） |
| 证据 | Siamycin I/II 是**三环肽（tricyclic）融合抑制剂**；耐药选择实验定位到 **gp160 氨基酸改变**（gp41 融合域），证实靶为包膜糖蛋白 |
| 关键结构 | Siamycin II 溶液结构（NMR，PMID 7787424）；无 siamycin–gp41 共结晶结构 |
| 文献 | Siamycin I 表征 PMID 8787894；Siamycin II NMR DOI 10.1007/BF00211754 / PMID 7787424；gp160 耐药突变 PMC163071 |
| 评级 | **低（同源迁移）**：无人类同源物，不能作「肽–人靶」正样本；仅可作「**肽–病毒蛋白**」预训练（先验学习肽–蛋白界面，再迁移到人靶），或直接**排除** |

---

## 六、汇总表

| 迁移对 | 细菌/病毒靶标 (UniProt) | 人同源物 (UniProt) | 结合界面保守性证据 | 关键 PDB | 文献 (PMID/DOI) | 评级 |
|---|---|---|---|---|---|---|
| RNAP 次级通道 | E. coli RpoC/β′ (**P0A8T7**) | POLR2A/RPB1 (**P24928**) | 桥螺旋+次级通道跨三域保守；amanitin 位点在漏斗/桥螺旋下方，与 MccJ25 相邻 | 6N60, 6N61, 1K83, 3CQZ | 30626643 / 10.1073/pnas.1817352116；11854475；19060880；20949106 | **高** |
| ClpC1–AAA+ | Mtb ClpC1 (**P9WPC3**) | CLPB (**Q9H078**)、CLPX (**O76031**) | AAA+ 核心中等保守；但 Lassomycin 结合 NTD 在人不保守 | 8IBO, 8YD1, 9IF4, 7XC5, 7US2, 7TTR | 24684906 / 10.1016/j.chembiol.2014.01.014；8IBO 相关 Int J Biol Macromol 2023 | **中** |
| 脯氨酰寡肽酶 | 细菌 POP (S9) | PREP (**P48147**) | S9 催化域保守；界面级复合物稀缺 | 1QFM, 1H2W（人 PREP） | Propeptin: J. Antibiot. (Kimura 等)；Biosci. Biotechnol. Biochem. 61:1754 | **低** |
| 跨拓扑环肽（CsA） | —（直接人靶） | PPIA (**P62937**) | 直接环肽–人蛋白复合物 | 1CWA, 2RMA, 1CYH | CsA–CypA 结构（Nature 等） | **高**（界面迁移） |
| 跨拓扑环肽（FK/Rap） | —（直接人靶） | FKBP12 (**P62942**) | 直接大环–人蛋白复合物 | 1FKF, 1FKD, 1FKJ, 1FAP, 1FKB | FKBP12–FK506/Rapamycin 结构 | **高**（界面迁移） |
| 跨拓扑环肽（其他） | —（直接人靶） | C3 / trypsin / TfR1 / NKp46 / ACE2 | 环肽/双环肽–人蛋白复合物 | 2QKI, 4GUX, 9GH7, 9H8R, 8BFW | 见各 PDB 条目 | **高**（界面迁移） |
| 病毒靶标 | HIV-1 Env gp160 (**P04578**) | 无同源物 | 无 | Siamycin II NMR（无共晶） | 8787894；7787424；PMC163071 | **低**（仅预训练/排除） |

---

## 七、推荐迁移训练集组合

**优先级排序（按迁移价值）：**

1. **RNAP 次级通道对（核心）**：
   - 正样本：MccJ25–RNAP（6N60）、capistruin–RNAP（6N61）→ 迁移到 α-amanitin–人 POLR2A（1K83 / 3CQZ）界面。
   - 理由：桥螺旋 + 次级通道跨域保守，且 amanitin 是**环肽–人 POLR2A 直接复合物**，构成「lasso–细菌RNAP ↔ 环肽–人RNAP II」的闭环迁移证据。

2. **跨拓扑环肽–人靶复合物（主扩充池）**：
   - CsA–PPIA（1CWA/2RMA）、FK506/雷帕霉素–FKBP12（1FKF/1FAP）、compstatin–C3（2QKI）、MCoTI-II–trypsin（4GUX）、双环肽–TfR1/NKp46/ACE2（9GH7/9H8R/8BFW）。
   - 用作「**环肽–人类蛋白**」界面正样本，附**拓扑标签**（lasso/cyclic/bicyclic/disulfide）做 domain adaptation。

3. **ClpC1–AAA+ 对（辅助，谨慎）**：
   - 仅迁移「肽–AAA+ 核心」的粗粒度特征（CLPB/CLPX），**不**迁移 Lassomycin 的 NTD 界面；可引入 ClpC1-NTD–Lassomycin（8IBO）作为「肽–细菌 AAA+」预训练。

4. **病毒靶标（可选，建议排除）**：
   - Siamycin–HIV Env 无人类同源物，不作「肽–人靶」正样本；最多用于「肽–蛋白界面」的通用预训练。

**最终推荐组合**：`{RNAP 对} ∪ {跨拓扑环肽–人靶复合物} ∪ {ClpC1 粗粒度 AAA+ 特征}`，**排除病毒靶标**，并在特征中显式编码肽拓扑类别以处理 lasso/环肽差异。

---

## 八、参考文献（PMID / DOI 索引）

1. Braffman et al., *Structural mechanism of transcription inhibition by lasso peptides microcin J25 and capistruin*, PNAS 116(4):1273–1278, 2019. **PMID 30626643**；DOI **10.1073/pnas.1817352116**（PDB 6N60/6N61）
2. Bushnell et al., *Structural basis of transcription: α-amanitin–RNA polymerase II cocrystal at 2.8 Å*, PNAS 2002（PDB 1K83）。**PMID 11854475**
3. Brueckner & Cramer, *Structural basis of transcription inhibition by α-amanitin*, PNAS 2008（PDB 3CQZ）。**PMID 19060880**
4. Weinzierl 等，*Conformational coupling, bridge helix dynamics and active site dehydration in catalysis by RNA polymerase*, 2010。**PMID 20949106**（桥螺旋附近：MccJ25/streptolydigin/α-amanitin 并列）
5. Gavrish et al., *Lassomycin, a ribosomally synthesized cyclic peptide, kills Mycobacterium tuberculosis by targeting ClpC1P1P2*, Chem. Biol. 21(4):509–518, 2014。**PMID 24684906**；DOI **10.1016/j.chembiol.2014.01.014**
6. ClpC1 NTD–Lassomycin 复合物晶体结构（Int. J. Biol. Macromol. 2023，PDB **8IBO**）
7. ClpC1P1P2 cryo-EM（PDB **8YD1** / **9IF4**）；人 CLPB ANK 域（PDB **7XC5**）；Skd3/CLPB AAA+（PDB **7US2/7TTR**）
8. CsA–CypA（PDB **1CWA/2RMA/1CYH**）；FKBP12–FK506（**1FKF/1FKD/1FKJ**）；FKBP12–雷帕霉素（**1FAP/1FKB**）
9. Compstatin–C3c（PDB **2QKI**）；MCoTI-II–trypsin（PDB **4GUX**）；双环肽–TfR1（**9GH7**）、–NKp46（**9H8R**）、–ACE2（**8BFW**）
10. Propeptin 发现（Kimura 等, J. Antibiot.）；Prolyl endopeptidase inhibitors from actinomycetes（Biosci. Biotechnol. Biochem. 61:1754, 1997）
11. Siamycin I 表征 **PMID 8787894**；Siamycin II NMR **PMID 7787424**（DOI 10.1007/BF00211754）；gp160 耐药突变 **PMC163071**

---
*注：UniProt / PDB 编号均以 web 检索命中为准；个别文献的年份/卷期以 PubMed 页面为权威来源。*

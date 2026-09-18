"""Write the target pre-survey document from the web search results."""
import os

ROOT = r'D:\deepseek_harness\prp49'

doc = r'''# 潜在成药靶点预调研（最终 10 对的靶）

> 2026-09-19 凌晨 · **网页快速检索**，目的是看这些靶**在文献里是否已有肽类/大环肽配体先例**
> ⚠️ **性质说明**：这是**广度优先的预调研**，依据是检索结果的标题与摘要级信息，
> **未逐篇精读全文**。任何一条都应在正式使用前回到原文核实。

---

## 一、为什么先查"有没有肽类配体先例"

本项目的候选来自 Lasso 肽平台。一个靶若**已有肽类/大环肽配体**，意味着：

1. **该靶的表面可被肽类占据** —— 不是只有小分子深口袋可及
2. **可能有共晶/冷冻电镜结构** —— 可用于对接与结合模式分析
3. **成药性判断更有依据** —— 肽类配体已证明可成药（如 MK-0616 已进入临床）

反之，若某靶从未有肽类配体报道，我们的候选就更像"模型的假设"而非"可类比的设计"。

---

## 二、逐靶调研结果

### 🔴 最重要的发现：EDNRB 的 Lasso 肽结构已发表

> **Structure of a lasso peptide bound ET_B receptor provides insights into the mechanism
> of GPCR inverse agonism** — *Nature Communications* **2025**, 16:3446
> DOI [10.1038/s41467-025-57960-x](https://www.nature.com/articles/s41467-025-57960-x) ·
> PMID [40263271](https://pubmed.ncbi.nlm.nih.gov/40263271/) · 冷冻电镜 **EMD-62275** · Nureki 组

- **这是 Lasso 肽结合人 GPCR 的结构工作**，机制是 **GPCR 反向激动（inverse agonism）**
- **EDNRB 正是本项目的靶之一**，RES-701-1/3 是我们的 **grade A/B 参考对**
- **项目意义（重要）**：
  1. **Lasso 肽靶向人蛋白有结构级别的先例** —— 不是只有细菌靶（如 MccJ25→RNAP）
  2. 该结构很可能就是 **9KDF** 的来源（我们的 grade-A 对、用于对接阳性对照）
  3. "反向激动"提示 Lasso 肽不只被动结合，**可调节受体功能** —— 对成药性叙事有利

### 最终 10 对的靶（逐个）

| 靶 | 是否有肽类/大环肽先例 | 检索到的线索 | 对项目的含义 |
|---|---|---|---|
| **JAK2** | 🟡 有环肽 mimetic | 基于 SOCS3 蛋白的**局部约束二甲苯环状 mimetic**（*Eur J Med Chem* 2025）| 该激酶界面可被环肽模拟，支持我们的候选 |
| **BCL2L1**（BCL-xL）| 🟢 **有环肽 + 结构** | PDB **7Y8D / 7Y99**：环肽 **cp1** 结合 BCL-xL | **有共晶可做对接验证**，是本批中结构证据最强者之一 |
| **ADORA2A** | 🟡 以**小分子**为主 | 腺苷/咖啡因类经典小分子配体；未检索到明确的肽类抑制剂 | 若我们的候选成立，属**较新颖**的方向（风险与机会并存）|
| **AURKA** | 🟢 **有靶向肽** | "Structure-based screening … novel **Aurora-A-targeting peptide** with antiproliferative activity"（*J Enzyme Inhib Med Chem*）| 已有肽类先例且带抗增殖活性 |
| **VHL** | 🟡 **小分子配体成熟** | VHL E3 连接酶小分子配体已广泛用于 **PROTAC**；未见肽类为主 | VHL 是"配体成熟"的靶，但**肽路线不主流** |
| **FURIN** | 🟡 拟肽抑制剂 | "Potent inhibitors of furin … containing **decarboxylated P1 arginine mimetics**"（*J Med Chem*）| 已有拟肽类抑制剂，**P1 精氨酸模拟**是关键化学 |
| **SRC** | 🟢🟢 **肽类大环抑制剂 + 机制研究** | "**Structural and Biochemical Basis for Intracellular Kinase Inhibition by Src-specific Peptidic Macrocycles**"（*Cell Chem Biol*）| **与本项目最匹配的靶之一**：胞内激酶 + 肽类大环 + 结构机制 |
| **ESR1**（ERα）| 🟢 **桥连肽 + 共晶** | "co-crystal structure of **peptide 12a** in complex with the **LBD of ERα**" | **有共晶结构**，支持肽类占据 LBD |
| **EGFR** | 🟢 **大环肽筛选** | "In Vitro Selection of **Macrocyclic α/β 3-Peptides** against Human EGFR"；EGFRp4 peptides | 大环肽筛选技术已用于该靶 |
| **CDK2** | 🟢 **肽抑制剂** | "Design of a **novel class of peptide inhibitors** of CDK/cyclin activation"（*JBC* 2005）| 肽类抑制 CDK2 有先例 |

**汇总**：10 个靶中 **6 个有明确的肽类/大环肽先例**（BCL2L1、AURKA、SRC、ESR1、EGFR、CDK2），
2 个以拟肽/环肽 mimetic 形式存在（JAK2、FURIN），2 个目前以**小分子**为主（ADORA2A、VHL）。

### 平台层面的旁证

- **大环肽已进入临床**：Merck 的 **MK-0616**（含氟-Trp 的大环肽，用于高胆固醇血症）
- **领域综述**："Nature-inspired macrocyclic peptides: Discovery and molecular engineering
  for drug development" —— 大环肽作为药物模态正在系统性推进

---

## 三、对项目的直接启示

### 3.1 支持我们候选合理性的证据

1. **SRC × Lassomycin**（我们 top 10 中**唯一 low-bias + 双信号一致**的对）
   → **恰好 SRC 是检索中"肽类大环抑制剂 + 结构机制"证据最强的靶之一** ✅ **独立佐证**
2. **BCL2L1 × Chaxapeptin** → BCL-xL 有环肽共晶（7Y8D/7Y99），**可立即做对接交叉验证**
3. **ESR1 / EGFR / CDK2 / AURKA** → 均有肽类先例，属"可类比设计"的靶

### 3.2 需要谨慎的靶

- **ADORA2A**：以腺苷类小分子为主，肽类先例不明 → 若候选成立是**新发现**，但失败风险也高
- **VHL**：小分子 PROTAC 配体成熟，**肽路线不主流** → 我们的候选在成药路径上不利

### 3.3 下一步可做的（低成本、高信息量）

| 动作 | 为什么值 |
|---|---|
| 用 **7Y8D/7Y99**（BCL-xL 环肽共晶）做对接阳性对照 | 检验我们的对接协议在"肽-人靶"上是否给出一致结论 |
| 精读 **SRC peptidic macrocycles** 论文 | 明确其结合位点与序列特征，与 Lassomycin 比对 |
| 核实 **Nat Commun 2025** 的 ET_B 结构与 **9KDF** 是否为同一数据 | 若是，我们的 grade-A 参考对有了正式文献出处 |
| 查 ESR1 的 **peptide 12a** 共晶（结构 ID 与结合模式）| 为"肽类可占据核受体 LBD"提供模板 |

---

## 四、本次预调研的局限（必读）

1. **仅依据检索结果的标题/摘要**，未逐篇精读；个别判断（如"未见肽类先例"）**不能当作否定证据**，
   只能说明"本轮检索未发现"
2. **检索于 2026-09-19**，文献进展快，结论有时效性
3. **未区分结合位点**：同一靶的不同位点成药性差别很大（如 JAK2 的 ATP 口袋 vs 变构位点），
   本轮只回答了"该靶是否有肽类配体"，**没有回答"与我们候选的结合位点是否相同"**
4. **部分链接为出版商站点**（Nature/EuropePMC），本轮抓取受 403 限制，未能读取全文
5. 结构化数据（PDB/EMDB ID）应作为**下一步核实**的入口，而非本轮的结论

---

## 五、引用（本轮检索命中的关键来源）

- Lasso 肽 × ET_B 受体结构：[*Nat Commun* 2025;16:3446](https://www.nature.com/articles/s41467-025-57960-x) ·
  [PubMed 40263271](https://pubmed.ncbi.nlm.nih.gov/40263271/) · [EMDB EMD-62275](https://www.ebi.ac.uk/emdb/EMD-62275)
- Src 肽类大环：[*Cell Chem Biol* — Structural and Biochemical Basis for Intracellular Kinase Inhibition by Src-specific Peptidic Macrocycles](https://www.sciencedirect.com/science/article/pii/S2451945616302549)
- BCL-xL 环肽结构：[PDB 7Y8D](https://pdbjlc1.pdbj.org/mine/summary/7Y8D)
- AURKA 靶向肽：[*J Enzyme Inhib Med Chem*](https://www.tandfonline.com/doi/pdf/10.1080/14756366.2026.2700842)
- EGFR 大环肽：[In Vitro Selection of Macrocyclic α/β 3-Peptides against Human EGFR](https://m2.mtmt.hu/api/publication/33166553)
- CDK2 肽抑制剂：[*JBC* 2005;280:13793](https://pubmed.ncbi.nlm.nih.gov/15649889/)
- FURIN 拟肽：[Potent Inhibitors of Furin … Decarboxylated P1 Arginine Mimetics](https://pubs-acs-org-443.webvpna.lzu.edu.cn/doi/pdf/10.1021/jm9012455)
- VHL 配体与 PROTAC：[Chem Soc Rev 2022](https://pubs.rsc.org/kr/content/articlehtml/2022/cs/d2cs00387b)
- JAK2 环状 mimetic：[*Eur J Med Chem* 2025 — Locally constrained xylene-based cyclic mimetics of SOCS3 protein](https://www.sciencedirect.com/science/article/pii/S0223523425009754)
- 大环肽综述：[Nature-inspired macrocyclic peptides](https://www.sciencedirect.com/science/article/pii/S1367593126000876)
'''
p = os.path.join(ROOT, 'docs', 'PRP49_潜在成药靶_预调研.md')
open(p, 'w', encoding='utf-8').write(doc)
print(f'wrote {p} ({len(doc.splitlines())} lines, {len(doc)} chars)')

import shutil
for rel in ['v2_2026-09-18_full', 'FINAL_2026-09-18_v2_full']:
    shutil.copy2(p, os.path.join(ROOT, 'releases', rel, 'PRP49_潜在成药靶_预调研.md'))
print('copied into both releases')

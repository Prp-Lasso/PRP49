# PRP49 论文插图集（10 张，待审核）

> 图 1–7 为 SVG（矢量，可无损缩放、可直接编辑）；图 8–10 为 PNG（matplotlib，200 dpi）。
> 每张图均给出**论文中的使用位置**。请审核后确认或提出修改意见。

## 图 1 · Lasso 肽拓扑：异肽键大环 + 穿线尾部

对比线性肽（易被蛋白酶降解）、普通环肽（刚性但无穿线）与 Lasso 肽（尾部穿过大环，形成 [1]rotaxane 拓扑，蛋白酶无法接近主链）。用于**绪论 1.1 与 1.2**。

![图 1 · Lasso 肽拓扑：异肽键大环 + 穿线尾部](figures/fig01_lasso_topology.svg)

*矢量源文件：`figures/fig01_lasso_topology.svg`*

---

## 图 2 · Lasso 肽的生物合成

前体肽 → lasso cyclase 催化形成内酰胺环 → C 端尾部穿环 → 成熟 Lasso 肽。底部标注核心矛盾：30 年仅约 50 条实验结构 vs 基因组挖掘出的数千条序列。用于**绪论 1.1**。

![图 2 · Lasso 肽的生物合成](figures/fig02_biosynthesis.svg)

*矢量源文件：`figures/fig02_biosynthesis.svg`*

---

## 图 3 · 已表征 Lasso 肽的靶点谱（细菌 vs 人源）

按证据等级着色（A = 共晶/冷冻电镜，B = 结合实验）。上半为细菌靶点（RNAP、ClpC1、核糖体），下半为人源靶点（ETB 受体、NPR1、GCGR）。用于**绪论 1.2**，支撑"Lasso 肽已可作用于人源靶点"。

![图 3 · 已表征 Lasso 肽的靶点谱（细菌 vs 人源）](figures/fig03_target_spectrum.png)

*位图：`figures/fig03_target_spectrum.png`*

---

## 图 4 · 通用结构预测工具的失效与 LassoPred 的解决

左：AlphaFold2/3、ESMFold 因套索结折叠与异肽键而失败；右：LassoPred 的两段式流程（分类器标注环/环区/尾部 + 构造器组装 3D），把结构覆盖从 <50 条扩展到 4,749 条。用于**绪论 1.3**。

![图 4 · 通用结构预测工具的失效与 LassoPred 的解决](figures/fig04_structure_prediction.svg)

*矢量源文件：`figures/fig04_structure_prediction.svg`*

---

## 图 5 · 本项目技术路线总览

八步流程：LassoPred 输出 → 分子对接 → Lasso/线性分类器 → 专用嵌入器 → 双线性注意力模型 → 样本稀缺处理 → 约束回归 → 筛选与开源。用于**方法总览**。

![图 5 · 本项目技术路线总览](figures/fig05_pipeline.svg)

*矢量源文件：`figures/fig05_pipeline.svg`*

---

## 图 6 · 双编码器双线性注意力架构

靶序列经 ESM-2（35M）、肽序列经 LassoESM（650M）编码，由 BAN 融合后**并联**输出分类头（是否结合）与回归头（结合强度）——两个头是并列的，不是串联。用于**方法 3.1**。

![图 6 · 双编码器双线性注意力架构](figures/fig06_architecture.svg)

*矢量源文件：`figures/fig06_architecture.svg`*

---

## 图 7 · 正样本稀缺与家族同源扩增

左：正样本从 352 扩增到 759；右：扩增后在更严的家族分组 CV 下**并未提升**（路线 A 仅 0.7481，低于基线 0.8142）。用于**方法 3.2 + 讨论（否定结果）**。

![图 7 · 正样本稀缺与家族同源扩增](figures/fig07_family_augmentation.png)

*位图：`figures/fig07_family_augmentation.png`*

---

## 图 8 · 四档验证协议阶梯（核心方法学图）

同一模型族在五种协议下的表现：随机初始化 0.594 → 分层 0.8145（有泄漏）→ 按肽分组 0.8333 → 按家族分组 0.8142（最严）→ 按靶分组 0.8738（同类分布）。误差棒为折间标准差。用于**结果 4.1**。

![图 8 · 四档验证协议阶梯（核心方法学图）](figures/fig08_cv_ladder.png)

*位图：`figures/fig08_cv_ladder.png`*

---

## 图 9 · 候选筛选漏斗 650 → 43 → 10

标注每级排除的数量与理由（已知靶 143 / 高偏置 200 / 无对接分 104 / 靶内排名>3）。底部列出全部硬性筛选规则。用于**结果 4.3**。

![图 9 · 候选筛选漏斗 650 → 43 → 10](figures/fig09_funnel.png)

*位图：`figures/fig09_funnel.png`*

---

## 图 10 · 最终 10 对候选及其融合权重稳健性

按融合得分排序，颜色表示稳健性（橙 = 全部 5 种权重下均入选、蓝 = 4/5、红 = ≤3/5），括号内为入选权重数。用于**结果 4.3**。

![图 10 · 最终 10 对候选及其融合权重稳健性](figures/fig10_final_candidates.png)

*位图：`figures/fig10_final_candidates.png`*

---

## 审核要点（请重点看这几处）

1. **图 1** 的穿线表达是否清楚——尾部是否看得出"穿过"大环（白色断口表示尾部在环前方）
2. **图 3** 的靶点清单是否要增删（当前 6 个细菌靶 + 4 个人源靶）
3. **图 5** 的八步流程是否与实际工作次序一致（尤其"对接"与"分类器"的先后）
4. **图 8** 的协议命名是否与论文正文一致（random init / stratified / peptide / family / target）
5. 配色与字号是否满足投稿要求（当前为英文标注，中文说明在正文图注中）

## 说明：一处事实更正

LassoPred 的合作方是 **Vanderbilt University（范德堡大学）**，位于**美国田纳西州纳什维尔**，
**不在德国**。论文作者与单位为：赵一雷（上海交大）、Zhongyue J. Yang 与 Ouyang Xingyu 等（Vanderbilt）、
A. James Link（普林斯顿）。绪论中已按此表述。

另：**LassoESM 是 UIUC + Vanderbilt（Mitchell 组）的工作**（Mi, Barrett, Mitchell, Shukla,
*Nat Commun* 2025），本项目是**采纳**该嵌入器，而非参与开发——绪论中已按此区分。
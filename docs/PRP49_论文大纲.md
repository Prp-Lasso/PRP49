# PRP49 结项论文 · 大纲

> 大纲正文约 1000 字 · 所有事实性陈述均已核实并附来源（见文末）
> 图片清单见第二节，共 10 张，已生成于 `docs/figures/`

---

## 拟定标题

**中文**：面向药物潜力挖掘的 Lasso 肽–蛋白相互作用预测：从专用嵌入器到候选短名单

**英文**：*Lasso peptide–protein interaction prediction for drug-potential mining: from a
domain-specific embedder to a prioritised candidate shortlist*

---

## 一、正文大纲

### 摘要（约 250 字）

Lasso 肽是微生物产生的一类含**异肽键**的特殊环肽，其大环与穿线尾部构成的**套索拓扑**
（[1]rotaxane / slipknot）赋予其显著的蛋白酶抗性与热稳定性，并已在抗菌、酶抑制等领域
展现出医药价值。然而该拓扑也使主流结构预测工具失效，并因实验数据稀缺而制约其效用预测。
本项目在指导老师团队前期工作 **LassoPred** 的基础上，接续开展**药物潜力导向**的 Lasso 肽
挖掘：先以分子对接建立结合可行性判据，构建 **Lasso 肽与普通链型肽的分类器**，
并引入专用嵌入器 **LassoESM**；进而搭建**双线性注意力（BAN）**神经网络。针对已报道正样本
稀缺导致的过拟合，采用**家族同源扩增**策略扩充训练样本；为实现有效初筛，建立**数值预测
（回归）模块**并系统比较多种**约束模式**，最终得到可用的预测模型与候选短名单，并已开源发布。

### 1 绪论

**1.1 Lasso 肽：一类含异肽键的特殊环肽**

Lasso 肽属于**核糖体合成与翻译后修饰肽（RiPPs）**。其核心特征有三：① N 端 α-氨基与
环内 Asp/Glu 侧链羧基形成**异肽键（isopeptide bond）**，闭合成大环（通常 7–9 残基）；
② C 端尾部**穿过**该大环，形成**套索（lasso）拓扑**；③ 部分成员含二硫键进一步锁定构象。
这种**穿线轮烷**结构在自然界罕见，30 年间仅约 **50 条**被实验解析，而生物信息学预测的
序列已达**数千条**。

**1.2 套索拓扑带来的医药相关特性**

穿线结构使 Lasso 肽兼具**大环肽的构象刚性**与**线性肽的序列可编码性**，并带来三项
医药相关优势：**蛋白酶抗性**（尾部被环包埋，主要攻击面仅剩暴露的 C 端）、
**热稳定性**、以及**对蛋白–蛋白界面的高亲和力**。已报道的活性实例横跨细菌与人源靶点：

| 肽 | 靶点 | 活性 | 证据等级 |
|---|---|---|---|
| MccJ25 | 细菌 RNAP β′ 亚基 | 抑制转录 | X 射线共晶（6N60）|
| Capistruin | 细菌 RNAP β′ | 抑制转录 | X 射线共晶（6N61）|
| Lassomycin | 结核分枝杆菌 ClpC1 | 抑制蛋白酶体 | X 射线共晶（8IBO/8IBP）|
| Lariocidin | 30S 核糖体 | 抗多重耐药菌 | 2025 年报道 |
| **RES-701-1/3** | **人 ETB 受体** | **拮抗** | **冷冻电镜（9KDF）** |
| **Anantin** | **人 NPR1** | **拮抗** | **放射性配体结合（Kd 0.6 µM）** |

**肠道微生物来源**是其中一个重要方向：MccJ25 通过**铁载体受体 FhuA** 被主动摄取进入
菌体（*Nat Chem Biol* 结构基础），Microcin Y 则被报道可调节肠道益生菌代谢并清除肠道
沙门氏菌。产业界亦已跟进——**Lassogen** 公司以"Lasso 肽疗法"为定位，管线覆盖
**感染、炎症与自身免疫疾病**，印证了该分子类别的成药潜力。

**1.3 特殊结构对结构预测的挑战**

Lasso 拓扑对计算预测构成双重障碍：不规则**套索结（lariat knot）**折叠，以及**非标准的
异肽键**。其结果是 **AlphaFold2、AlphaFold3 与 ESMFold 均无法准确预测** Lasso 肽结构。
本项目指导老师**赵一雷**（上海交大）与**美国范德堡大学（Vanderbilt University）**
Zhongyue J. Yang 团队及**普林斯顿大学** A. James Link 团队合作开发的 **LassoPred**
（*Nat Commun* 2025）解决了这一问题：以分类器标注环/环区/尾部，再由构造器组装 3D 结构，
将已知结构覆盖从不足 50 条扩展至 **4,749 条**计算模型，并显著优于上述通用工具。
> ⚠️ 需注意：范德堡大学位于**美国田纳西州纳什维尔**（非德国）。

**1.4 效用预测的挑战与本项目定位**

结构的解决并不等于功能的解决。通用蛋白质语言模型在 Lasso 肽相关任务上表现不佳，
为此出现了专用嵌入器 **LassoESM**（Mi 等，*Nat Commun* 2025），其在底物兼容性与
RNAP 抑制活性预测上均有提升。**本项目正是在这一基础上，把目标从"结构预测"推进到
"药物潜力挖掘"**：不仅判断一条 Lasso 肽是否折叠正确，而是预测它**能否与特定人类靶点结合**。

### 2 前驱工作

**2.1 分子对接**：以 AutoDock Vina 建立 Lasso 肽 × 靶点的对接流程，处理大环肽的
盒尺寸与刚性约束问题，建立结合可行性判据。

**2.2 Lasso 肽 vs 普通链型肽的分类器**：构建序列层面的判别模型，用于从候选肽库中
区分真正的 Lasso 拓扑（本项目中作为"拓扑合理性"过滤器）。

**2.3 专用嵌入器采纳**：引入 LassoESM 作为肽编码器，与靶标编码器（ESM-2）配对。

### 3 方法

**3.1 双编码器与双线性注意力（BAN）**：靶序列经 ESM-2、肽序列经 LassoESM 编码，
经 BAN 双向注意力融合后接分类头（是否结合）与回归头（结合强度）。

**3.2 正样本稀缺与家族同源扩增**：公开报道的 Lasso 肽–人靶正样本极少，直接训练
严重过拟合。**在不违背初筛样本要求**的前提下，以**家族同源**方式扩充正样本，
并按**家族分组**做交叉验证以确保评估无泄漏。

**3.3 数值预测模块与约束模式**：为支撑有效初筛，搭建回归模块预测结合强度。
针对亲和力数据跨度大、易发散的问题，系统比较多种**输出约束模式**，确定可用配置。

### 4 结果与讨论

**4.1 分类性能**：报告四档递进验证协议下的性能（分层 / 按肽分组 / 按家族分组 / 按靶分组），
以量化"数据泄漏"对指标的贡献。

**4.2 回归与约束**：报告各约束模式的 RMSE 与相关性，说明"收敛但未达实用精度"的边界。

**4.3 候选短名单与开源发布**：经拓扑过滤、成药性靶点筛选与双信号（模型 + 对接）融合，
产出候选短名单并**开源发布**，同时给出**稳健性分级**。

### 5 结论

建立 Lasso 肽–人蛋白互作的预测流程，界定其适用边界（分布内筛选可行、分布外发现不可行），
并交付可复现的模型、数据与候选清单。

---

## 二、插图清单（10 张，均已生成）

| # | 图名 | 类型 | 文件 |
|---|---|---|---|
| 1 | Lasso 肽拓扑：异肽键 + 大环 + 穿线尾部（对比线性肽与普通环肽）| 结构示意 | `fig01_lasso_topology.svg` |
| 2 | Lasso 肽生物合成与穿线过程 | 流程示意 | `fig02_biosynthesis.svg` |
| 3 | 已表征 Lasso 肽的靶点谱（细菌靶 vs 人源靶）| 分类图 | `fig03_target_spectrum.svg` |
| 4 | 通用结构预测工具的失效与 LassoPred 的解决 | 对比图 | `fig04_structure_prediction.svg` |
| 5 | 本项目技术路线总览 | 流程图 | `fig05_pipeline.svg` |
| 6 | 双编码器 BAN 模型架构 | 架构图 | `fig06_architecture.svg` |
| 7 | 正样本稀缺与家族同源扩增 | 数据图 | `fig07_family_augmentation.svg` |
| 8 | 四档 CV 协议阶梯 | 柱状图 | `fig08_cv_ladder.png` |
| 9 | 候选筛选漏斗：650 → 43 → 10 | 漏斗图 | `fig09_funnel.png` |
| 10 | 最终 10 对候选及其稳健性分级 | 条形图 | `fig10_final_candidates.png` |

---

## 三、需与组员确认的事项

1. 联合署名与章节分工（本组员工作见 `PRP49_与组员工作整合方案.md`）
2. 组员的 lpiLasso 回归模型在**分组 CV** 下的复评结果（见 `给组员的复评说明_分组CV.md`）
3. 组员的 `translasso`/`TailBackbone` 结构生成是否纳入本文（可与 LassoPred 做对照）

---

## 四、参考文献（已核实）

1. Ouyang X, Ran X, Xu H, Al-Abssi R, **Zhao YL**, Link AJ, Yang ZJ. *LassoPred: a tool to
   predict the 3D structure of lasso peptides.* **Nat Commun** 16:5497 (2025).
   doi:10.1038/s41467-025-60544-4
2. Mi X, Barrett SE, Mitchell DA, Shukla D. *LassoESM: a tailored language model for enhanced
   lasso peptide property prediction.* **Nat Commun** 16:8545 (2025). doi:10.1038/s41467-025-63412-3
3. Ouyang X, Ran X, Zhu D, Yang ZJ. *Predicting lasso peptide structure with LassoPred.*
   **Methods Enzymol** 730:91-106 (2026). doi:10.1016/bs.mie.2026.01.037
4. *Structural basis for hijacking siderophore receptors by antimicrobial lasso peptides.*
   **Nat Chem Biol** (PMC3992131)
5. *Lasso Peptides — A New Weapon Against Superbugs* (review). PMID 40943111
6. *Lasso peptides: A focus on therapeutic index.* PMID 40289244
7. *Systematic mining of the human microbiome identifies antimicrobial peptides with diverse
   activity spectra.* PMID 37973865
8. *Microcin Y ... effectively clear gut Salmonella.* Int J Biol Macromol (2024)
9. *The lasso structure, biosynthesis, bioactivities and potential applications of Microcin J25.*
   (2023) — 综述
10. Lassogen — Lasso Peptide Based Therapeutics（在研适应症：感染 / 炎症 / 自身免疫）
11. 本团队数据来源：PDB 6N60 / 6N61 / 6N62（RNAP×MccJ25/Capistruin）、9KDF（ETB×RES-701-3）、
    8IBO / 8IBP（ClpC1×Lassomycin）、4CU4（FhuA×MccJ25）

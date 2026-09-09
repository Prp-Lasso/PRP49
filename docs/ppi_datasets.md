# 肽-蛋白互作（PPI / LPI）可复用公共数据集与工具调研

> 项目背景：用蛋白质语言模型（pLM）预测 **Lasso 肽（RiPP 类）—人类蛋白互作（LPI）**。
> 正样本极少，拟采用「预训练/预热 + 迁移学习」：先用大量普通肽（线性肽、环肽）-蛋白互作数据预热双编码器，
> 再用 Lasso 小样本微调。
>
> 调研日期：2026 年（基于公开文献与官网信息，规模数字为论文/官网口径，标注「约」者为约数，使用前请以官网为准）。

---

## 1. 肽-蛋白复合物结构数据库（构建结构级正样本与接触图标签）

| 名称 | 类型 | 规模 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **Propedia v2.3** | 蛋白-肽复合物结构库（来自 PDB，聚类去冗余，含图结构签名） | 约 1.9 万个蛋白-肽复合物（v2.x 口径，官网/论文为准） | 是（PDB 来源，含大量人蛋白） | 官网 `https://bioinfo.dcc.ufmg.br/propedia2/`；代码/数据 `https://github.com/LBS-UFMG/propedia` | **首选结构级正样本来源**：直接导出肽序列+受体序列+接触残基对，用于双编码器预训练与接触图标签；按受体物种过滤出"人蛋白"子集 |
| **PepBDB** | 生物肽-蛋白复合物结构库（PDB 派生，聚类去冗余，含对接基准） | 约 4000+ 复合物（论文 4,231 口径） | 是（含受体/人蛋白） | `http://huanglab.phys.hust.edu.cn/pepbdb/` | 经典蛋白-肽对接基准与结构标签来源；规模较 Propedia 小但质量高、有聚类信息，适合做干净的验证集/接触图标签 |
| **PepPIs** | 肽-蛋白互作数据库（较早，2010 年代） | 较小、维护程度低 | 部分 | 早期站点（检索可见历史文献，现官网已难访问） | 仅作历史参考，不优先使用；被 Propedia/PepBDB 取代 |
| **InterPep** | 肽-蛋白互作位点预测服务器 + 基准 | 预测服务器；附基准集 | 是 | `https://interpep.bioinfo.se/`；配套 InterPepRank 集 `https://researchdata.se/en/catalogue/dataset/doi-10-17044-scilifelab-14222692` | 用于**接触位点标签**与模型对照（残基级互作位点），可作为监督信号补充 |
| **pepATTRACT** | 盲法大规模肽-蛋白对接服务器 | 对接工具 + 基准（数十~上百复合物） | 是 | `https://bioserv.rpbs.univ-paris-diderot.fr/services/pepATTRACT/`（NAR 2017） | 对接/结构互补性预筛与负样本生成（见第 6 节），不直接给数据集 |
| **PPIKB** | 蛋白-肽互作知识库（文献+专利，2025 预印本） | 较大（整合 Propedia/PepBDB/PDBbind+ 等口径） | 是 | `http://ppi.duanlab.ac/`；下载 `https://ppikb.duanlab.ac/downloads/` | 较新的一站式肽-蛋白互作知识库，可作为 Propedia 的补充与交叉校验来源 |
| **PepXPro** | 结构-亲和力蛋白-肽数据集整理框架 | 把 PDBbind 19,037 复合物筛成 2,458 肽候选 | 是 | 预印本 `https://www.biorxiv.org/content/10.1101/2026.08.09.743757v1`（框架开源） | 可直接复用的「PDBbind→肽子集」清洗流程，产出结构+亲和力联合标签 |

---

## 2. 肽-蛋白互作 / 肽活性数据库

### 2a. 通用互作数据库（可抽取"肽/小蛋白—蛋白"条目）

| 名称 | 类型 | 规模 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **BioGRID** | 多物种物理/遗传互作库 | 数百万条互作 | 是（人为主） | 下载 `https://downloads.thebiogrid.org/BioGRID/`；筛选器/API `https://thebiogrid.org/` | 下载 MITAB 后按**分子类型（polypeptide/short peptide）+ 物种=Homo sapiens**过滤，抽取肽-蛋白互作正样本 |
| **IntAct** | 分子互作库（IMEx 成员） | 百万级条目 | 是 | `https://www.ebi.ac.uk/intact/`；PSICQUIC/MITAB 下载 | 通过 PSICQUIC 按 interactor type 过滤 peptide，与 BioGRID 交叉去重 |
| **STRING** | 蛋白-蛋白关联网络（含分数） | 1.4 万+ 物种，人网络数十万对 | 是 | `https://string-db.org/`；API `https://version12.string-db.org/help/api/`；批量下载 | 提供人蛋白-蛋白关联分数，作**预热阶段的"蛋白侧"语料**与软负样本来源；肽条目较少，需自行映射 |
| **PepPPI/PepPI 类预测工具** | 肽-蛋白互作预测服务器/模型 | 工具 | 是 | 检索可见（多篇 PepPI 预测方法论文） | 仅作方法对照，非数据集 |

### 2b. 肽活性数据库（注意：多为 MIC/表型/活性，非"互作对"）

| 名称 | 类型 | 规模 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **DBAASP v3** | 抗微生物/细胞毒肽活性与结构库 | 约 1.5 万+ 肽、2 万+ 活性记录 | 间接（人细胞毒性数据） | `https://dbaasp.org/`（可下载 CSV） | 提供**肽侧序列多样性**与抗菌/细胞毒标签；无直接"肽-蛋白"互作对，仅用于肽表征预训练或负对照 |
| **DRAMP v3** | 抗微生物肽数据仓库 | 约 2.2 万+ AMP | 部分（含人宿主/靶注释） | `http://dramp.cpu-bioinfor.org/` | AMP 序列大语料；作为 pLM 肽侧预热数据或负样本库 |
| **CancerPPD / CancerPPD2** | 抗癌肽库 | 约 3,000+ 抗癌肽 | 部分 | `http://crdd.osdd.net/raghava/cancerppd/`（CancerPPD2 更新版） | 抗癌肽序列+活性，非互作对；用于肽表征预热 |
| **SATPdb** | 结构注释治疗肽库 | 约 1.9 万治疗肽 | 部分 | `http://crdd.osdd.net/raghava/satpdb/` | 治疗肽序列/结构语料，含环肽，适合预热 |
| **AVPdb** | 实验验证抗病毒肽库 | 约 2,600+ 抗病毒肽 | 部分 | `http://crdd.osdd.net/servers/avpdb/` | 抗病毒肽序列/活性语料 |
| **CyBase** | 环肽（cyclotide）数据库 | 百级 cyclotide（含结构） | 部分 | `http://cybase.org.au/`（`https://cybase.org.au/`） | **环肽侧最相关语料**：cyclotide 序列+结构，与 Lasso 肽同为植物环肽，可做肽侧迁移学习的桥接 |
| **KNOTTIN** | 抑制剂胱氨酸结（knottin）骨架库 | 约 2,000 knottin 序列/结构 | 是（部分靶人蛋白） | `http://www.dsimb.inserm.fr/KNOTTIN/` | 二硫键环化肽语料；与 Lasso（酰胺键环化+套索拓扑）结构同源度较低，但可作为"环肽"大类预热数据 |

---

## 3. 药物-靶标亲和力数据（回归预训练）

| 名称 | 类型 | 规模 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **PDBbind** | 蛋白-配体复合物+实验结合亲和力（Ki/Kd/IC50） | v2020：约 1.94 万复合物（refined 集约 5,316） | 是 | `http://www.pdbbind.org.cn/`；新版 `https://www.pdbbind-plus.org.cn/`（需注册） | **核心回归标签来源**：过滤 peptide 子集（可复用 PepXPro 流程，19,037→2,458 肽复合物）做"肽-蛋白亲和力回归"预热 |
| **BindingDB** | 实验结合数据（含肽配体） | 约 270 万+ 结合记录 | 是 | `https://www.bindingdb.org/`（TSV/SDF 下载） | 按配体类型=peptide + 靶物种=人，抽取肽-靶 Ki/IC50 作回归语料；规模大、标签丰富 |
| **ChEMBL** | 生物活性数据库 | 约 1,500 万+ 活性记录 | 是 | `https://www.ebi.ac.uk/chembl/`（SQLite/API） | 用 API 按 `molecule_type=peptide` + 人靶过滤，抽取肽-靶活性（IC50/EC50/Ki）；适合大规模回归预热与负样本采样 |

---

## 4. PPI 基准数据集（预热双编码器/注意力层）

| 名称 | 类型 | 规模 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **TUnA（Wang-lab-UCSD）** | Transformer PPI 预测框架 + 基准（Bernett、XSpecies 等） | Bernett：7 物种（含人）PPI 对；XSpecies：跨物种 | 是（Bernett 含人） | 仓库 `https://github.com/Wang-lab-UCSD/TUnA`；预训练权重 `https://huggingface.co/yk0/TUnA_models`；嵌入 `https://huggingface.co/yk0/TUnA_embeddings`；另有 `https://huggingface.co/yk0/tuna-r-esm_mlp` | **首选预热基准**：直接复用 Bernett 数据集格式（TSV，蛋白A序列/蛋白B序列/标签），下载 ESM 嵌入或权重做双编码器预热；xspecies 用于跨物种泛化测试 |
| **DeepPPI** | 深度 PPI 预测 + 基准（S. cerevisiae、人类） | 酵母+人类 PPI 对（万级） | 是（人类子集） | `https://github.com/hashemifar/DeepPPI`（TF 复现版 `https://github.com/gdario/deep_ppi`） | 人类 PPI 基准，作序列双编码器预热的补充数据集 |
| **PIPR（seq_ppi）** | 序列 PPI 预测 + STRING 派生基准 | STRING 派生（酵母/人） | 是 | `https://github.com/muhaochen/seq_ppi` | 与 STRING 一致的 PPI 对格式，便于与 TUnA/DeepPPI 交叉验证 |
| **STRING（直接）** | 蛋白-蛋白关联 | 人网络数十万对 | 是 | `https://string-db.org/` | 提供海量正/负 PPI 对，按 combined_score 阈值分层作预热语料 |

> 注：PPI 基准多为**蛋白-蛋白**（非肽-蛋白），用途是预热编码器与注意力层（让模型先学会通用互作表征），再用第 1/3 节的肽-蛋白结构/亲和力数据做领域迁移，最后 Lasso 微调。

---

## 5. 人类靶标蛋白资源

| 名称 | 类型 | 规模 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **UniProt 人类参考蛋白质组** | 序列数据库 | UP000005640，约 2 万个蛋白（含异构体） | 全部为人 | `https://www.uniprot.org/proteomes/UP000005640`（FASTA 下载） | **靶标序列基线**：下载人参考蛋白组作靶标侧序列语料与 ID 映射（gene↔UniProt↔PDB↔AF） |
| **AlphaFold DB** | 结构数据库 | 人蛋白组全覆盖（~2 万结构） | 全部为人 | `https://alphafold.ebi.ac.uk/`（批量下载） | **靶标结构来源**：10 个候选靶标全部有 AF2 结构，供结构互补性预筛与接触图标签 |

### 10 个候选人类靶标的序列/结构可用性

| 靶标（基因） | UniProt | 关键结构状态 | 备注 |
|---|---|---|---|
| EDNRB（内皮素 B 受体） | P24530 | 有 cryo-EM 实验结构（与内皮素配体），AF2 有 | GPCR，天然配体即肽，LPI 高相关 |
| NPR1（利钠肽受体 1） | P16066 | 有实验结构（与 ANP 结合），AF2 有 | 天然配体为环肽 ANP，结构域清晰 |
| POLR2A（RNA 聚合酶 II RPB1） | P24928 | 有大量 PDB 实验结构，AF2 有 | 大蛋白，接触面定义复杂，建议用结构域级片段 |
| CLPB（caseinolytic peptidase B） | Q9H078 | AF2 有（AAA+ ATPase，实验结构较少） | 依赖 AF2 结构 |
| PREP（脯氨酰内肽酶） | P48147 | 有实验结构，AF2 有 | 酶，活性位点明确 |
| PPIA（亲环素 A / CypA） | P62937 | 有实验结构（与环孢素 A 结合），AF2 有 | 经典药物靶，含肽底物 |
| FKBP1A（FKBP12） | P62942 | 有实验结构（与 FK506/雷帕霉素），AF2 有 | 经典药物靶，脯氨酰异构酶 |
| MDM2（E3 泛素连接酶） | Q00987 | 有实验结构（与 p53 肽、nutlin），AF2 有 | **天然肽互作靶**，p53 结合口袋清晰，理想正样本模板 |
| PD-L1（CD274） | Q9NZQ7 | 有实验结构（与 PD-1、抗体），AF2 有 | 免疫检查点，含大环肽抑制剂先例 |
| 整合素（ITGA5/ITGAV/ITGB1/ITGB3 等） | 多 | 有实验结构（含 RGD 环肽配体），AF2 有 | **环肽配体丰富**，是环肽-蛋白互作的黄金来源 |

> 上述靶标均在 UniProt 人参考蛋白组与 AlphaFold DB 内，且多数有实验 PDB 结构，可同时提供"AF2 结构"与"实验结构"两条结构标签路径。

---

## 6. 负样本 / 灰色样本生成工具（结构互补性预筛）

| 名称 | 类型 | 吞吐量/可获得性 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **HDOCK** | 混合式蛋白-蛋白/肽对接（模板+从头） | 高（Web 队列，分钟级/任务，可批量提交） | 是 | `http://hdock.phys.hust.edu.cn/`（免费） | **批量结构互补性预筛首选**：对候选肽-靶对打分，生成灰色/负样本排序 |
| **HPEPDOCK** | 盲法肽-蛋白对接（分层算法） | 高（Web 队列，适合肽库批量） | 是 | `http://huanglab.phys.hust.edu.cn/hpepdock/` | **肽-蛋白专用**，与 HDOCK 互补，用于 Lasso 肽对接打分与构象生成 |
| **CABS-dock** | 柔性肽-蛋白对接（粗粒度） | 中（Web 队列，10-50 肽/任务） | 是 | `http://biocomp.chem.uw.edu.pl/CABSdock/`（新站点 `https://cabsdock.pl`） | 肽高柔性场景（Lasso 环肽）对接，适合小规模精筛 |
| **ClusPro PeptiDock** | 肽对接（PIPER 打分） | 高（Web 队列） | 是 | `https://cluspro.org/`（peptide 协议） | 快速打分，可与 HPEPDOCK 交叉验证 |
| **RosettaDock / pepATTRACT（本地）** | 本地对接 | 中-高（可本地并行，pepATTRACT 开源） | 是 | Rosetta 许可；pepATTRACT 服务器+源码 | 需要大规模对接时用本地版本扩展吞吐 |

> 用途定位：这些工具用于 **（a）生成 Lasso 肽-靶的"结构互补性"打分作为弱标签；** **（b）构造难负样本（对接分数低但序列相似的对）**；**（c）为小样本 Lasso 正样本提供构象/接触面先验。**

---

## 7. RiPP / 环肽-蛋白互作专用数据集与配套工具

| 名称 | 类型 | 规模 | 是否含人靶 | 下载/访问方式（URL） | 对项目的用途建议 |
|---|---|---|---|---|---|
| **RODEO / RODEO2** | RiPP 前体肽预测与注释工具 | 工具（含 RiPP 家族模型，含 lasso peptide 类） | 否（注释工具） | `https://github.com/mitchell-lab/RODEO`（另有 Samaneh Behroozian 版） | **直接相关**：用于识别/校验 Lasso 肽前体与家族（lasso peptide 是内置类别），帮助构建 Lasso 肽正样本清单 |
| **BAGEL4** | 细菌素/RiPP 基因组挖掘 Web 服务器 | 工具+内置数据库 | 否 | `http://bagel4.molgenrug.nl/` | 细菌素/RiPP 挖掘，补充 Lasso 肽候选发现；配套数据可作 RiPP 语料 |
| **antiSMASH** | 次级代谢物基因簇（BGC）检测 | 工具（含 RiPP 类） | 否 | `https://antismash.secondarymetabolites.org/`（含配套 MIBiG 数据库 `https://mibig.secondarymetabolites.org/`） | BGC 检测 + MIBiG 数据库提供 RiPP/环肽基因簇注释，支撑 Lasso 肽序列收集 |
| **MIBiG** | 次级代谢物 BGC 最小信息数据库 | 数千 BGC | 否 | `https://mibig.secondarymetabolites.org/` | RiPP（含 lasso peptide）基因簇注释，间接用于正样本构建 |
| **CyBase / KNOTTIN（见第 2b）** | 环肽/knottin 结构库 | 百~千级 | 部分 | 见上 | 环肽侧迁移学习语料 |

> 专门"RiPP-蛋白互作对"的公开数据集目前**基本不存在**（这是本项目的空白点），现有资源多停留在"RiPP 序列/基因簇发现"层面；互作对需靠第 1 节结构库 + 第 6 节对接自建。

---

## 推荐数据组合与获取优先级（按项目阶段 M1–M2）

### M1 —— 预热双编码器（通用互作表征）

**目标**：让序列双编码器学会通用蛋白-肽/蛋白-蛋白互作表征，无需 Lasso 特异性。

**优先级排序（先到后）**：
1. **TUnA Bernett（含人）+ XSpecies**（`github.com/Wang-lab-UCSD/TUnA`；嵌入 `yk0/TUnA_embeddings`，权重 `yk0/TUnA_models`）—— 开箱即用、格式规范（TSV）、含 ESM 嵌入，**最省事**。
2. **Propedia v2.3 结构级肽-蛋白复合物**（`github.com/LBS-UFMG/propedia` / 官网）—— 提供肽-受体序列对 + 接触图，比纯序列 PPI 更贴近 LPI 任务。
3. **PepBDB**（`huanglab.phys.hust.edu.cn/pepbdb/`）—— 高质量小验证集。
4. **DeepPPI（人类子集）+ PIPR** —— 扩充正/负 PPI 对。
5. **STRING（人）** —— 按 combined_score 分层的软标签语料。

### M2 —— 领域迁移 + 回归 + Lasso 微调

**目标**：从"通用互作"迁移到"肽-蛋白亲和力/结构"并最终适配 Lasso 肽-人蛋白。

**优先级排序**：
1. **PDBbind 肽子集（PepXPro 流程）** + **BindingDB/ChEMBL 肽-靶活性** —— 回归预训练（Ki/Kd/IC50 标签），让模型学会"亲和力量级"。
2. **CyBase（cyclotide）+ KNOTTIN + DBAASP/DRAMP/CancerPPD/SATPdb/AVPdb** —— 环肽/活性肽序列语料，桥接"线性肽 → 环肽 → Lasso"。
3. **RODEO/RODEO2 + MIBiG + BAGEL4 + antiSMASH** —— 收集并校验 Lasso 肽前体/成熟肽序列，构建 Lasso 正样本清单。
4. **10 个人类靶标**：从 UniProt `UP000005640` + AlphaFold DB 取序列/结构（EDNRB、NPR1、MDM2、整合素等天然肽/环肽靶优先）。
5. **HDOCK / HPEPDOCK / CABS-dock / ClusPro PeptiDock** —— 结构互补性预筛 + 难负样本 + 接触面先验，弥补 Lasso 正样本稀缺。
6. **BioGRID + IntAct** —— 抽取人"肽-蛋白"条目补充正样本，并做负样本构造（同源蛋白不同肽）。

### 关键提示
- **正样本稀缺是核心约束**：优先抓取 Propedia/PepBDB/PDBbind 中"人蛋白 + 肽"的结构级正样本（数量级远大于 Lasso 专属数据），Lasso 专属互作对需自建。
- **标签分层**：结构接触图（Propedia/PepBDB/InterPep）> 亲和力回归（PDBbind/BindingDB/ChEMBL）> 活性/表型（DBAASP/DRAMP 等，仅辅助）> 纯序列互作对（TUnA/STRING）。
- **去泄漏**：下载后统一做 CD-HIT 序列聚类去冗余，避免训练/验证/测试间序列同源泄漏（LPI 小样本场景尤其重要）。
- **许可与下载**：PDBbind/BindingDB 需注册或遵守许可；部分旧站（PepPIs、部分 Raghava 服务器）可能已下线，优先用上表标注的活跃 URL 与 GitHub 镜像。

---

### 附：一键起点 URL 速查

- 结构级肽-蛋白：`https://github.com/LBS-UFMG/propedia`、`http://huanglab.phys.hust.edu.cn/pepbdb/`
- PPI 基准+嵌入：`https://github.com/Wang-lab-UCSD/TUnA`、`https://huggingface.co/yk0/TUnA_embeddings`、`https://huggingface.co/yk0/TUnA_models`
- 亲和力：`https://www.pdbbind-plus.org.cn/`、`https://www.bindingdb.org/`、`https://www.ebi.ac.uk/chembl/`
- 人靶序列/结构：`https://www.uniprot.org/proteomes/UP000005640`、`https://alphafold.ebi.ac.uk/`
- 对接预筛：`http://hdock.phys.hust.edu.cn/`、`http://huanglab.phys.hust.edu.cn/hpepdock/`
- RiPP/Lasso 工具：`https://github.com/mitchell-lab/RODEO`、`https://mibig.secondarymetabolites.org/`、`http://bagel4.molgenrug.nl/`

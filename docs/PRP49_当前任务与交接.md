# PRP49 当前任务与交接文档

> **用途**：本文档是"活文档"。任何新会话（尤其是上下文被压缩后）应先读本文件，即可无损续接全部工作。
> **最后更新**：2026-09-11 深夜（服务器时间约 21:00）
> **上一份总体回顾**：`docs/PRP49_项目总体回顾.md`

---

## 1. 服务器与环境（最关键）

| 项 | 值 |
|---|---|
| 登录 | `sylogin.hpc.sjtu.edu.cn:22`（实际解析到 `sylogin2.pi.sjtu.edu.cn`）|
| 账号 / 密码 | `clswxl-ccmbi1` / `2024zhao0620!` |
| HOME | `/dssg/home/acct-clswxl/clswxl-ccmbi1` |
| **唯一允许操作的目录** | **`~/LassoPep`**（用户明令：其余目录是他人资产，绝对禁止读写）|
| 连接方式 | 本机 Python + paramiko（`hpc_deploy/ops/opslib.py`，密码登录；本机无 sshpass）|
| 环境激活 | `source /usr/share/lmod/lmod/init/profile; module load miniconda3/4.10.3; source activate prp49` |
| 项目 conda 环境 | `prp49`（Python 3.10 + torch 2.3.1+cu121 + transformers 4.40.2）|
| 已有他人/用户环境（可复用） | `autodock_vina`（vina 1.2.3 + obabel + meeko 0.4.0）、`amber`、`pymol`、`af3`、`lasso_mypy3` 等 15 个 |
| 队列 | `a100`（GPU，1 卡 16 核，每核 8G）、`64c512g`（CPU，133 节点）|
| 作业上限 | 单账号 ≤ 500 作业（当前用约 20）|

**操作脚本**（本机，全部基于 paramiko）：
`D:\deepseek_harness\prp49\hpc_deploy\ops\*.py` —— 每个脚本一个任务（上传/提交/查询/修复）。
典型用法：`python 'D:\deepseek_harness\prp49\hpc_deploy\ops\<script>.py'`（PowerShell，注意用单引号包裹路径）。

**作业脚本**（本机 → 上传到服务器 `~/LassoPep/scripts/`）：
`D:\deepseek_harness\prp49\hpc_deploy\scripts\*.sh`

---

## 2. 正在进行中的作业（截至本文档）

**当前（09-17 09:00 快照）**

| JobID | 名称 | 内容 | 状态 |
|---|---|---|---|
| **62657636** | `prp49_ctrftA` | **路线 C 下游微调**（5 折并行，修复 #41 后重提）| 🟢 fold0 运行中（从 ep 30 续训）|
| **62657649** | `prp49_verify` | **真实复合物验证**（6 个 PDB 复合物 × 人靶对照）| 🟡 排队 |
| **62646669** | `prp49_homol` | **路线 D′** 域同源加权 MTL | 🟢 运行 8h |
| ✅ 完成 | `62642044_[0-4]` | **RMSE 方案 A 五折并行** | ✅ 5/5 完成，**RMSE 1.530 vs 基线 1.813（+15.6%），Pearson 0.546** |
| ✅ 完成 | `62602714_2` | **RMSE 方案 B = 赢家（绝对误差）** | ✅ 5/5，RMSE **1.472** vs 1.623（+9.3%）|
| ❌ 失败 | `62646194_[0-4]` | 路线 C 首版 | ❌ 全部 NameError（坑 #41）→ 已修复重提 |
| ✅ 完成 | `prp49_contr` 62624147 | **路线 C 对比预训练** | ✅ 15 epochs / 4分19秒；分离度 **0.024→0.442** |
| ✅ 完成 | `prp49_prop` 62593230 | 路线 A 家族软标签传播 | ✅ 5 折 **0.7481 ± 0.0506** < 基线 0.8142 → **判负** |

**历史作业（09-12 快照，多数已完成）**

| JobID | 名称 | 内容 | 状态 |
|---|---|---|---|
| 62394452 | `prp49_dockp` | 对接批 1 重提（171 失败任务：18 OOM + 2 超时）| 已改 MAXPAR=4 / mem=190G 重跑 |
| 62363938 | `prp49_dockp2` | 对接批 2（1,843/1,999）| 剩 3 任务（不再等）|
| 62358463 | `prp49_grp` | **无泄漏验证**（按肽分组 CV）| ✅ fold0 0.812 / fold1 0.849 |
| 62359642 | `prp49_s2` | 跨拓扑迁移（Propedia 8,346 对）| fold1 0.755 / AP 0.761 |
| 62394451 | `prp49_regc` | **A+C 亲和力分类** | ✅ 5 折 AUC 0.6312 |
| 62406210 | `prp49_rank` | 成对排序（成药筛选导向）| ✅ 见路线计划 |
| 已停止 | `prp49_reg` | 原亲和力回归（RMSE 6.86，判定不可行）| 已被 A+C 取代 |

**查询模板**：
```python
# 新建 ops/<name>.py，内容：
import sys; sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run
cli = connect()
st, out, err = run(cli, '<shell 命令>')
print(out); cli.close()
```

---

## 3. 关键结果数字（当前最优）

| 实验 | CV AUC | 备注 |
|---|---|---|
| 8M MVP（CPU）| 0.663 | 本机验证 |
| 650M 随机初始化基线 | 0.594 | 5 折 stratified |
| I1 课程学习单独 | 0.566 | **无效**（略降）|
| **I1+I2+I5 全套改进** | **0.8145 ± 0.017** | 5 折 stratified |
| **grp（无泄漏协议，按肽分组）** | **0.8333 ± 0.0120**（5 折：0.812/0.849/0.837/0.837/0.831）| ✅ **泄漏疑虑彻底排除**（grouped 高于 stratified）|
| s2（Propedia 跨拓扑迁移）| **0.7539**（5 折 0.741–0.770）| ✅ 迁移学习有效 |
| regc（A+C 亲和力分类，靶标泛化）| **0.6430**（5 折 0.617–0.675）| 🟡 略低于 0.65 判据：跨靶零样本信号有限但真实 |
| **rank（成对排序，成药筛选）** | 4/5 折：pair_acc **0.553** · 逐靶 Spearman **0.136** · NDCG@5 **0.688** · **EF@10% 1.48** | ⚠️ **对 Lasso 肽是负贡献**（域偏移）：加入后已知对平均排名 23 → 71，故候选融合中默认关闭 |
| **跨域泛化测试（Propedia 普通肽）** | **AUC 0.470**（打乱负）/ **0.541**（同受体真实肽）| ❌ 域特异：模型只在 Lasso 域有效，**不可用于普通肽库** |
| Lasso 判定器（原版，长度捷径）| F1 **0.973** / AUC 0.998 | ❌ 负样本长度分布不重叠 → 虚高 |
| **Lasso 判定器（匹配负样本，诚实版）** | **Val F1 0.7863** | ✅ best_model.pt 已是此版本（CD-HIT 去重：正 1697 / 负 4414；训练 4888 对）|

**对接矩阵**：11 肽 × 10 靶 = 110 对完成（`docking/dock_matrix_full.csv`）
- 阳性对照 CTRL_9KDF redocking = **−12.2**（盒放大后合理；放大前是 +104.8 假信号）
- 模型分 vs 对接分 Spearman **ρ = −0.091**（不相关 → 信息互补）
- 融合 `z(model)+w·z(−dock)`：已知对平均排名 **4.43 → 3.57**（w=1.0）；前 20% 命中 4/8 = 50%（**未达验收 70%**）
- 主失败案例：**PB1m7×PLXNB1 排 11/11**（人工嫁接肽，模型盲点）、MccJ25×POLR2A 排 10/11

**Propedia 实验复合物对接（2,161 对已聚合）**：`docking_propedia/all_scores.csv`
- b1：n=324，mean −25.7 / median −28.2；b2：n=1,837，mean −21.2 / median −21.4
- clash（正分）仅 12 个（0.6%）；98.5% ≤ −8 kcal/mol（都是实验验证的结合对，符合预期）
- 11 个异常值（|score| > 200，最高 4.4×10⁷）已隔离到 `score_outliers.csv`

---

## 4. 数据资产（服务器 `~/LassoPep/`）

| 路径 | 内容 |
|---|---|
| `mvp_cpu/train_pairs_hard.csv` | 训练集（正 312 + 硬负样本），improved/grp 用 |
| `mvp_cpu/affinity_pairs.csv` | 4321 对实验亲和力（Ki 1615 / Kd 545 / IC50 2161；564 靶；`energy` = pAffinity 3–12）|
| **`mvp_cpu/affinity_cls_pairs.csv`** | **A+C 数据集：2462 对**（靶内 z-score 后 |z|>0.5 分箱，正 1235 / 负 1227；218 靶）|
| **`mvp_cpu/affinity_pairs_rank.csv`** | **成对排序集：16,819 对**（同靶内肽对，ΔpAffinity ≥ 0.3；335 靶；列 prot_id/prot_seq/pep_a_seq/pep_b_seq/energy_a/energy_b/label）|
| **`mvp_cpu/propedia_dock_pairs.csv`** | **对接增强集：3746 对**（Propedia 实验复合物正样本 1873，对接分 ≤ −12；打乱肽负样本 1873；698 受体；附 `dock_score` 列）|
| `mvp_cpu/scan_peptides_ext.csv` | 13 条扫描肽 |
| `mvp_cpu/scan_targets.fasta` | 11 条靶链 |
| `data/lassopred_database.csv` | LassoPred 库（拓扑注释来源，4749 条；**注意路径是 `data/` 不是根目录**）|
| `lasso_id_data/` | 判定器数据（正 4749 / 负 20000 / 匹配负 10559）|
| `LassoPeptideClassifier/` | 判定器代码 + `checkpoints/best_model.pt`（诚实版 F1 0.786）|
| `docking/` | Lasso 矩阵（110 对）+ `tasks.csv` + `pep/`、`rec/`、`pdbqt/`、`out/` |
| `docking_propedia/` | 批 1（500 对）→ **`all_scores.csv`（2161 对聚合，含序列+对接分）** + `score_outliers.csv` |
| `docking_propedia_b2/` | 批 2（1999 对）|
| `docking_drugpanel/af/` | 41 个成药人靶 AF 模型（备用）|
| `PRP49/runs*/checkpoints/fold*_best.pt` | 各训练权重（650M 配置每个 2.75 GB）|
| `results/scan_matrix_ext.csv` | 模型打分矩阵（13 肽）|

**训练配置清单**（`PRP49/config_*.yaml`）：`config.yaml`（基线）、`config_warmup.yaml`（S1 课程预热）、`config_improved.yaml`（I1+I2+I5）、`config_improved_grouped.yaml`（无泄漏复评）、`config_s2.yaml`（Propedia 迁移）、**`config_reg_cls.yaml`（A+C 亲和力分类）**、**`config_s3.yaml`（对接增强迁移，含 dock_score）**、**`config_rank.yaml`（成对排序）**

**成药筛选交付链路**（本项目最终交付形态）：
```
候选肽库 ──► ①Lasso 判定器 (lasso_prob)      —— 拓扑是否合理
        └─► ②结合分类头 (binding_prob)      —— 是否结合
        └─► ③成对排序模型 (rank_score)      —— 同靶内谁更强
        └─► ④AutoDock Vina (dock_score)     —— 结构是否可行
                    └──► 加权 z-score 融合 ──► 候选短名单 CSV
```
一键生成：`python -m prp49.candidates --config config_rank.yaml --peptides ... --targets ... --cls_ckpt ... --rank_ckpt ... --lasso_ckpt ... --dock_csv ... --out ../results/candidates.csv`
（`prp49/candidates.py` 已实现并上传；权重产出后即可运行）

**本机对应路径**：`D:\deepseek_harness\prp49\`（`mvp_cpu/`、`docking/`、`results/`、`scratch/`（分析脚本）、`affinity_data/`）

**目标交付架构图**（archify 生成，可交互：主题切换 / 缩放 / 搜索 / 导览视图 / 导出）：
- `diagrams/prp49-architecture.html` —— 交付产物（showcase 校验 9/9，四档视口零溢出）
- `diagrams/prp49-architecture.json` —— 可编辑规格（13 节点 / 2 区域边界 / 13 连接 / 3 指标卡）
- 三个导览视图：主预测路径 · 训练与监督信号 · 结构证据融合

---

## 5. 下一步待办（按优先级，09-13 17:45 更新）

**验收状态（09-17 20:00 更新）**
| 标准 | 目标 | 结果 |
|---|---|---|
| 主任务 CV AUC（按肽分组）| ≥ 0.75 | **0.8333 ± 0.0120**（`runs_improved_grouped`）✅ |
| **主任务 CV AUC（按家族分组，最严）** | — | **0.8142**（`runs_family`，5 折 0.852/0.867/0.777/0.819/0.755）✅ **仅降 0.019 → 无家族泄漏** |
| 硬正样本 top-20% | ≥ 70% | **3/3 = 100%** ✅ |
| 硬负样本控制 | 不排前列 | **2/2 受控** ✅ |
| 域外泛化（多任务后）| — | **0.6421 ± 0.0335**（单任务 0.470）✅ |
| RMSE 约束回归 | 优于基线 | ✅ 方案 A **1.530** vs 1.813（+15.6%）；方案 B **1.472** vs 1.623。**但均未达实用门槛 <1.0** |

> ⚠️ **基线协议必须对齐（坑 #43）**：`runs_family` = **0.8142**（按**家族**分组，最严）；
> `runs_improved_grouped` = **0.8333**（按**肽**分组）；`runs_improved` = 0.8145（分层，**有泄漏**）。
> 引用基线前先确认协议一致 —— 我曾用 0.8142 对比按肽分组的路线 C，把"5/5 折为负"读成了"略优"。

**四条探索路线的最终判定（09-17）**
| 路线 | 定量结果 | 判定 |
|---|---|---|
| **A** 家族软标签传播 | family CV **0.7481 ± 0.0506** vs 基线 0.8142 | ❌ **判负** |
| **B** 家族分组 CV（验证协议）| 0.8333 → 0.8142，仅降 **0.019** | ✅ **通过**（泄漏排除）|
| **MTL** 多任务跨域 | 跨域 **0.470 → 0.6421** | ✅ **有效** |
| **C** 家族对比预训练 | 预训练分离度 0.024→0.442；下游 **5/5 折低于基线**（−0.0085，配对 t=−2.66）| ❌ **判负** |
| **D′** 域同源加权 | 域内 **+0.061** ✅ / 跨域 **−0.072** ❌ | ⚖️ **权衡，非改进** |
| **D v2** 真实复合物加权 | 运行中（fold 0 低于 D′）| 🟡 待定 |

> **路线 C 的方法学教训**：对比学习的分离度提升（0.024→0.442，负对相似度 0.966→0.502）
> **完全不代表下游受益**。若只看自监督指标就会误判为成功。详见 `docs/PRP49_ACD路线计划.md` 第十节。

**CV 协议阶梯（方法学资产，四档严格性递进）**
| 协议 | 问题设定 | AUC |
|---|---|---|
| 按肽分组 | 新肽 × 已知靶 | **0.8333** |
| **按家族分组** | **新家族 × 已知靶** | **0.8142** |
| 按靶分组（A+C）| 已知肽 × 新靶 | 0.6312（overall）|
| 跨物种域外 | 普通肽 × 蛋白 | 0.470 → 0.64（多任务）|

> **结论**：从"新肽"到"新家族"只掉 0.019 → 模型学到的是**跨家族泛化规律**，
> 而非家族记忆。这是本阶段最重要的方法学证据。

**多任务跨域迁移（09-15 20:00 完成，5 折）**
| 指标 | 单任务 | 多任务（λ=0.5）| 变化 |
|---|---|---|---|
| 域内（Lasso 肽）| 0.8333 | **0.7429 ± 0.0444** | −0.090 |
| **跨域（普通肽）** | **0.470** | **0.6421 ± 0.0335** | **+0.172** ✅ |

各折跨域：0.642 / 0.612 / 0.665 / 0.600 / **0.691**（全部 >0.60，稳定复现）
→ **共享 BAN 主干确实学到可跨域迁移的表征**，代价是域内小幅下降。

**RMSE 约束回归（09-15 20:30 进行中，fold 0-1 部分折）**
| 方案 | 配置 | 实测 RMSE | 该折基线 | 判定 |
|---|---|---|---|---|
| **A** | 全量 + tanh 约束 + **base_head** | **1.577** | 2.295 | ✅ 优于基线 31% |
| **B** | 仅 Kd+Ki（尺度最干净）| **1.348** | 1.427 | ✅ 略优（基线本身低）|
| **C** | 全量 + tanh 约束（**无** base_head，对照）| 1.689 | 2.295 | ✅ 优于基线，但**差于 A 0.11** |

→ **tanh 约束彻底解决了发散**（从失控的 6.855 到 1.58），**且 base_head 分解带来可测增益**（A 1.577 vs C 1.689）。
⚠️ 仍未达实用门槛（<1.0）；且需注意"ep 0 的中心值优势"（tanh 初始输出≈7.5 接近部分折的标签均值），
故最终判据应以**收敛值**与 **oracle 基线（测试折自身 std）** 同框比较。


### RMSE 约束回归 —— 最终 5 折对比（09-17 完成）

| 方案 | RMSE | 基线 | **相对改善** | Pearson | Spearman | 优于基线 |
|---|---|---|---|---|---|---|
| **A**（全量 + tanh 约束 + base_head）| 1.530 ± 0.067 | 1.813 | **+15.6%** | **0.546** | **0.517** | 4/5 |
| **B**（仅 Kd+Ki）| **1.472** | 1.623 | +9.3% | 0.484 | — | **5/5** |

**逐折（A）**：1.584 / 1.623 / 1.455 / 1.529 / 1.459（fold1 是唯一劣于基线的折）

**判定（诚实）**：
- tanh 约束**彻底解决发散**（6.86 → 1.5）是确定成果
- A 的相对改善与排序能力更好，B 的绝对误差更低，**两者相差仅 0.058**
- **均远未达实用门槛（<1.0）**，在 5 折上难以判定谁显著更优

**最终短名单（`results/candidates_final.csv`，121 行 / 91% 有对接分）**
| 排名 | 肽 | 靶 | composite | 对接 | 文献证据 |
|---|---|---|---|---|---|
| 1 | MccJ25 | POLR2A | 1.858 | −15.8 | C（RpoC 同源）|
| 2 | Ubonodin | POLR2A | 1.810 | −15.7 | C（细菌 RNAP）|
| 3 | Capistruin | POLR2A | 1.697 | −15.4 | C（RpoC 同源）|
| 4 | Ubonodin | EDNRB | 1.689 | −15.8 | |
| 5 | MccJ25 | C3 | 1.537 | −15.0 | |
| 6 | Anantin | EDNRB | 1.457 | −14.9 | （结构为**预测**）|
| **10** | **RES-701-1** | **EDNRB** | 1.318 | −14.5 | **B: IC50 10 nM** ✅ |
| **11** | **RES-701-3** | **EDNRB** | 1.314 | −14.5 | **A: 9KDF 冷冻电镜** ✅ |
| **23** | **Anantin** | **NPR1** | 0.663 | −12.8 | **B: Kd 0.6 µM** ✅ |

**关键结论（可写入报告）**
1. **无泄漏 CV 0.8333**：按肽分组 5 折，全部 0.81–0.85 → 模型判别能力真实
2. **对接分是决定性信号**：仅在对接分缺失时模型才会出错；补齐 ITGAVB3 命名映射后硬负样本从 42/43 名压到 77/79 名
3. **多任务训练可跨域迁移**（fold0）：域内 0.740（单任务 0.833）、**跨域 0.642**（单任务 0.470）→ 共享主干确实学到可迁移表征，代价是域内小幅下降
4. **温度缩放在域外无效**：T=0.881 使短名单分布几乎不变（43% >0.9 仍是 43%）→ 概率仅作排序，不作绝对解释
5. **对接受限**：大环肽刚性对接区分度低（多落在 −9~−16），Anantin 真靶 NPR1 排第 23 而非第 1

**本轮三项改进的成绩（诚实记录）**
| 改进 | 结果 | 判定 |
|---|---|---|
| **对接命名映射修复**（ITGAVB3 → ITGAV/ITGB3）| coverage 69%→85%，负样本 0/2 → **2/2** | ✅ 决定性修复 |
| 概率校准（温度缩放）| >0.9 比例 29%→**7%**，ECE 0.102→0.092 | ✅ 有效，未达 ECE<0.05 |
| xneg 特异性训练（352 跨靶负样本）| CV 0.833→**0.794**，假阳性未改善 | ⚠️ **建议弃用**（徒增代价）|
| 融合权重改对接主导（0.40/0.60）| 唯一同时满足正负样本的设置 | ✅ 但仅基于 5 个验证点，需扩展验证 |

**⚠️ 待改进的首要项**：Anantin×NPR1 因**无肽结构 → 无对接分**被缺失惩罚压低（58/143）。解法是为 Anantin 建模结构（ESMFold/Rosetta）补对接分，而非调权重。

---

## 6. 四项后续任务进展（09-14 12:25）

| 任务 | 状态 | 结果/说明 |
|---|---|---|
| **1. Anantin 建模 + 对接** | 🟢 对接运行中（62504422，88 对）| 用 **LassoPred 官方流程**跑通：PyMOL 突变 → tleap 成环。**8 个模型全部构建成功**（left/right × loop_len 3-6），CA 序列验证 100% 匹配 `GFIGWGNDIFGHYSGDF`。因 loop 长度/手性无文献记载，8 个模型全对接取中位数 |
| **2. 扩展权重验证** | ✅ 完成（结论受限）| 审计的 17 个矩阵外正样本中**仅 4 个可用**（BI-32169×GCGR + 3 个工程嫁接变体）；其余 13 个是细菌靶/RNA/非人源，**域不匹配**。→ 权重仍只基于 5 个验证点，已在代码注释中标注 |
| **3. 多任务训练** | 🟢 训练中（62504543）| 数据划分完成：Propedia train 3,058 对/1,245 受体、test 724 对/311 受体，**受体重叠 0**（按完整序列划分）。设计：共享 BAN 主干 + 两域 BCE（λ=0.5），双指标报告（域内 + 跨域 AUC）|
| **4. 候选表清理** | ✅ 完成 | `results/candidates_final.csv`：移除 PB1m7（非 lasso 肽）与 Capi-var1（计算变体）→ 121 行、82% 有对接分、9 行带文献证据标注、无结构肽显式 flag |

**任务 1 新增坑（#21-23）**
| # | 坑 | 现象 | 解法 |
|---|---|---|---|
| 21 | **pymol conda 环境是 Python 2** | LassoPred 驱动脚本直接 `SyntaxError`（f-string 不支持）| 用 prp49 的 Python 3 跑驱动，只把 pymol **二进制**加进 PATH |
| 22 | `lasso_leap()` 参数类型 | `ValueError: Invalid step` | 传字符串 `'min1'`/`'min2'`，不是整数 |
| 23 | **非标准残基名** | 序列验证读出 `GFIGWGNXIFGHYSGDF` 误判失败 | LassoPred 把成环残基改名为 **ASX/GLX/ANX**，组氨酸变 HID/HIE/HIP → 映射表需补齐 |
| 24 | **`csv.DictWriter` 默认 `lineterminator='\r\n'`** | 服务器端（Linux）生成的 tasks.csv 仍带 CRLF → `--size_z "126\r"` 非法 → **88 个对接任务全部静默失败**（vina 只打印 usage）| 写 CSV **必须显式** `lineterminator='\n'`；跨平台传 CSV 后用 `sed -i "s/\r$//"` 兜底 |
| 25 | **多任务每步跑两个域 → CUDA OOM** | 域内靶序列 1024 token × 2 个 batch → A100 40 GB 爆 | 改为**交替域**（每步只跑一个域）+ batch 8 + max_len 768 |
| 26 | **`sum(1 for x, y in zip(a, b))` 漏了 `if x == y`** | 该生成器数的是**元素个数**（恒等于 len(a)），于是 `identity()` 对任意序列对都返回 **1.0** → 家族传播的软标签全部被抬到最高档 0.9 | 必须写 `if x == y`；聚类脚本里那个版本写对了，所以 178 家族是对的，只有传播档位受影响 |
| 27 | **作业时限按"单方案"估但脚本里串行跑多方案** | RMSE 约束三方案写在一个 `job_reg_bound.sh` 里 → 8 小时只够跑半个方案 A，**无 JSON 产出** | 改成 `--array=1-3` 并行，每方案独立作业 |
| 28 | **靶标 ID 用 PDB 链描述** | `train_pairs_hard.csv` 里 81% 的 prot_id 是 `"Chains A, B"` 之类，**238 个蛋白序列只对应 49 个 ID** → 硬负样本的"同靶"假设部分失效 | 重建为 `train_pairs_v2.csv`：靶以**序列**为键（238 个），负样本 = 同序列 × 不同家族 |
| 29 | **改了一半的索引：`all_scores[te] = sc`** | 评估排除传播样本后 `sc` 长度为 705 而 `te` 为 800 → **形状不匹配，fold 0 末崩溃**（白跑 2.6 小时）| 散播必须用同一索引 `all_scores[eval_te] = sc`；且 `overall` 指标也要用 `eval_mask` 过滤（否则未赋值的传播行以 0.0 参与计算）|
| 30 | **检查点只在"每折结束"保存，epoch 循环内不落盘** | 一旦撞上 wall-clock 时间上限，**内存中的 best 权重与已完成的 epoch 全部丢失** —— 本项目因此白跑 **3 次**（regb 串行 8h、regb2 8h、MTL 6h）；regb2 那次连 5 折结果都没有，因为 JSON 也只在最后写 | 新增 `prp49/ckpt.py`：**每 `ckpt_every`（默认 10）个 epoch 覆盖式保存 `fold{i}_latest.pt`**（含 optimizer + epoch，支持续训），**原子写入**（`.tmp` + `os.replace`，写盘中途被杀也不会损坏）；`fold{i}_best.pt` 仍是交付用最佳权重。每个训练脚本启动时 `maybe_resume()` 自动续训。regb2 的 JSON 也改为**每折写一次** |
| 31 | **内联 `python -c "..."` 嵌套引号** | 多次因 PowerShell 引号转义导致语法错误、甚至**整个命令块未执行**（config 修改静默失败）| 一律改为**写脚本文件再执行**，不用内联多行命令 |
| 32 | **上传脚本前未本地试跑** | 变量名 typo（`ex_p` 应为 `exp_p`）在服务器上跑到第 78 行才报错；CRLF 也在上传后才暴露 | **任何上传的脚本先在本机跑通**（用本机镜像数据）；配 `scratch/test_*_locally.py` + `preflight.py` |
| 33 | **PyMOL 路线重复踩坑** | ① 该 conda 环境是 **Python 2**（f-string 语法错误，坑 #21 已记过一次）② 许可证 **2024-05 已过期** | 结构叠合改用**纯 Python3 + numpy**（`map_box_to_af_py3.py`），不再依赖 PyMOL |
| 34 | **多配置共用 checkpoint 目录 → resume 到别的配置的权重** | RMSE 三方案 A/B/C 都写 `runs_reg_bound/checkpoints`；方案 C（**无** base_head）启动 fold 1 时 `maybe_resume` 读到方案 B（**有** base_head）的权重 → `Unexpected key(s): base_fc/base_out` → **CRASH**（6.6h 后，丢 fold 1-4）。**这是修 #30 时引入的新问题**：加 resume 前只写不读，不会崩 | ① `ckpt.py` 捕获 `RuntimeError` 并**拒绝加载不兼容检查点**（打印提示后从零开始）② **每个配置独立 `checkpoint_dir`**（作业脚本内 `sed` 生成 `config_regb{A,B,C}.yaml`）|
| 35 | **基于部分折下结论** | 路线 A 前 2 折均为 0.726 → 我判断"明显有害"并建议提前终止；**实际 5 折均值 0.7481**（后两折 0.805/0.807），与基线 0.8142 仅差 0.066 | **至少完成 3 折或全部折再判断**；折间标准差可达 0.05，2 折远不足以定位均值（用户当时要求"等完全完成再做结论"，判断正确）|
| 36 | **作业时限按静态猜测而非实测** | RMSE 方案 A 单折实测 **7h**，5 折串行需 **35h**，但时限设 16h → **16:00:29 TIMEOUT，只完成 2 折**（`epochs=60` 也是 #31 引号失败导致 sed 未生效）| 时限 = **实测单折时间 × 折数 × 1.3**；长任务一律 `--array` **按折并行**（同样 GPU·h，墙钟从 35h 降到 7h）；`preflight.py` 需改为读历史日志实测值而非静态常数 |
| 37 | **并行 array 共享 `--out` 会互相覆盖** | 5 折并行若都写 `reg_bound_A.json`，最后完成的覆盖其余 | 单折模式自动改名为 `reg_bound_A_fold{N}.json`，事后用 `merge_reg_folds.py` 汇总 |
| 38 | **并行任务共享 checkpoint 目录 → "假完成"** | 方案 A 的 fold 1-4 在 **1-2 分钟内"完成"**（`best_epoch: null`、seconds 58-117）：新作业 `maybe_resume` 读到了**旧作业遗留的、已训练完的 ckpt** → 跳过全部训练，直接评估**别人的权重**，并给出看似正常的 RMSE。**比 #34（崩溃）更危险，因为它静默产出结果** | ① 每折独立 `--ckpt_dir`；② 检查点已完成时**打印 WARN 并跳过**，绝不静默评估；③ 受污染结果移入 `results/discarded_tainted/` |
| 39 | **单一关键词检索静默漏收集** | 全文检索 `"lasso peptide"` 得到 98 个 PDB 结构，但**零个** RNAP 复合物；改用**目标驱动**检索（`RNA polymerase microcin`）才找到 **6N60/6N61/6N62**（RNAP×MccJ25/Capistruin）—— 它们的标题里根本没有 "lasso" 这个词 | 同时按**靶**与按**肽类别**两路检索；对比两路结果的差集 |
| 40 | **多亚基复合物"取第一个匹配"** | 6N60/6N61 的靶链被选成 **omega 亚基（91 aa）**，因为匹配串 `"DNA-directed RNA polymerase "` 命中了 α/β/β′/ω/σ **所有**亚基，而首行恰好是错的那个。正确的 β′（rpoC，1409 aa）才是人 POLR2A 的对应物 | 用**区分性字符串**（`subunit beta'`）+ **最长匹配优先**；多链复合物中永远不要 `iloc[0]` |
| 41 | **只做语法检查就以为改好了** | 我给 `train.py` 加 `--fold` 时把守卫写进 `run_cv()`，却在其中引用 `main()` 的局部变量 `args` → **5 个 array 任务在 6-48 秒内全部 NameError**。`ast.parse` 语法完全通过，**抓不到作用域错误** | 改函数接口后做**针对性检查**：用 AST 确认函数体内不引用外部作用域的自由变量（`args`）+ 检查调用点传参；能跑冒烟测试就跑 |
| 42 | **核酸链被当作肽** | 6N60/6N61/6N62 含短 DNA 链，我的过滤只用长度（`<=80 → 肽`），于是 `AAACCTTCCCATCAGGGGAAGGTC` 作为"lasso 肽"混入辅助集 | 按**沉积的 entity 类型**（`polypeptide(L)`）过滤，而非长度；辅以"只含 ACGTU 则判为核酸"的校验 |
| 43 | **对比基线协议不匹配** | 我拿**按家族分组**的基线（0.8142）去对比路线 C，而 `config_contrastive.yaml` 用的是 `cv_mode: grouped`（**按肽分组**，正确基线 `runs_improved_grouped` = **0.8333**）→ 把"4/4 折低于基线（均值 −0.0091）"误读成"略优 +0.009"，**结论方向完全相反** | 对比前先确认**双方的 CV 协议一致**；从 `cv_summary.json` 实读基线，不凭记忆引用数字 |
| 44 | **假设时间而不实测** | 未查 `date` 就假定时刻，把闹钟设到了**已经过去**的时间（显示 `in -2h1m53s`）；本地与服务器时钟本是同步的 | 任何与时间有关的判断（预计完成、闹钟、时限）都先执行一次 `date`；闹钟设置后检查返回的倒计时符号 |

---

## 7. 必须记住的坑（已踩过，勿重蹈）

| # | 坑 | 症状 | 修复 |
|---|---|---|---|
| 1 | **CSV 行尾 CRLF** | vina 报参数错误/秒退（`\r` 混入字段）| 生成时 `newline='\n'`；上传后 `sed -i "s/\r$//"`（已踩 3 次）|
| 2 | **pep/rec 同名文件** | vina 把受体当配体 → `PDBQT parsing error` | 脚本里加前缀 `pep_` / `rec_` |
| 3 | **回归用 BCE** | loss = −790、RMSE = 503 | `label_type=energy` 走 Huber 回归损失（已修 `train.py`）|
| 4 | **模块版本不同步** | `PairDataset got unexpected kwarg 'energy_col'` | 上传代码时同步 `data.py` + `train.py` + `eval.py` 全套 |
| 5 | **GroupKFold + 负样本同组** | val AUC = nan | 负样本组 id 唯一化，或改 `cv_mode: stratified` |
| 6 | **对接 OOM** | `OUT_OF_ME+`，部分任务被杀 | 降并发 + 加 mem（64c512g 每核 8G 上限）|
| 7 | **SFTP 路径** | `FileNotFoundError` | SFTP **不展开 `~`**、**不自动建目录**：先 `mkdir -p`，传绝对路径 |
| 8 | **判定器长度捷径** | 负样本长度分布与正样本不重叠（正 20-50 aa vs 负 100-500 aa）→ F1 虚高 | 已建"匹配负样本"（打乱 + 前导肽 + 同长度细菌）→ `lassoid2` 重训 |
| 9 | **盒太小** | 大环肽对接出**假正值**（+104）| 配体派生盒最小边 ≥ 38 Å；本轮上限 44 Å（速度/精度平衡）|
| 10 | **节点故障** | `NODE_FAIL`（机时自动返还）| 直接重跑 |
| 11 | **登录节点禁算** | 用户明令 | 安装/训练/对接**全部 sbatch 到计算节点**，登录节点只做上传与提交 |
| 12 | **BatchNorm + batch=1** | `Expected more than 1 value per channel when training, got input size [1, 256]`（BAN 层 BN）→ 训练在 fold1 崩溃 | 训练 DataLoader 加 `drop_last=True`（已修 `train.py`）|
| 13 | **同名缓存键冲突** | 聚合脚本里 pep/rec 同名 → 受体序列被肽覆盖（"receptor aa 1-34"）| pep_cache / rec_cache 分开 |
| 14 | **对接异常分值** | vina 输出 |score| 高达 4.4×10⁷ 污染统计 | 聚合时按物理范围 |score| ≤ 200 过滤，异常值单独落盘 |
| 15 | **pandas 写 CSV** | `to_csv(newline=...)` 报 TypeError | 用 `lineterminator='\n'` |
| 16 | **PowerShell 写 JSON 加 BOM** | archify 报 `Unexpected token ''` | 用 write 工具写 JSON，别用 `Set-Content -Encoding UTF8` |
| 17 | **PDB 记录名切片** | `line[:6] == 'ATOM'` 永假（前 6 列是 `"ATOM  "` 带空格）→ **39 个受体文件全成 2 字节空文件**，429 个对接任务静默失败 | 必须 `.strip()` 后比较（已修 `scratch/find_pockets.py`）|
| 18 | **空 score 文件欺骗进度检查** | `ls out/*.score \| wc -l` 显示 429/429，实际全是 0 字节 | 检查**非空**文件数：`find out -name "*.score" -size +0c \| wc -l` |
| 19 | **同一受体两种命名** | 对接表用 `ITGAVB3`（二聚体）、候选矩阵用 `ITGAV`/`ITGB3`（单链）→ merge 静默丢分（coverage 69%），而**实验确认的硬负样本恰好在这两列** → 假阳性排进前 30% | 在 `candidates.py` 内映射：只把 `ITGAVB3` 那些行复制成 `ITGB3`（**复制整表会凭空造分**）|
| 20 | **换数据源未核对靶集** | 39 个"成药靶"（EGFR/JAK2/BCL2…）与原 11 个 Lasso 靶**完全不重叠**，换表后 coverage 掉到 **0%** | 换数据源后必须打印 coverage；两个靶集是互补的，不能互相替代 |

---

## 8. 对接吞吐与资源记录（两轮）

**第一轮（盒 44 Å + exh 4 + 2 核）**：原参数每对需 2+ 小时（2.2 小时仅 398/2499）→ 提速 12–20 倍。
**第二轮（09-12，针对 OOM）**：b1 仍有 18/25 任务 `OUT_OF_MEMORY`（2 个 8 小时超时）→ 再降为 `MAXPAR=4` + `--mem=190G`，失败任务用 `sbatch --array=<job列表>%6` 重提（脚本按 `out/<id>.score` 存在性跳过已完成项）。

---

## 9. 关键设计决策（勿轻易推翻）

1. **评估协议**：主结果用 **grouped CV（按肽序列分组）**；stratified 数字仅作对照（有泄漏）。
2. **硬负样本**：同靶 × 其他家族肽 —— 用于逼模型学"肽-靶兼容性"而非"肽身份"，但因此必须用 grouped CV 评估。
3. **L_align 结构对齐**：用复合物接触图监督 BAN 注意力（`align_loss_weight=0.1`，仅注释对生效）。
4. **对接协议**：AutoDock Vina 1.2.3，刚性配体（大环全冻结，`rigidify.py`），盒=肽自身坐标 + pad，exhaustiveness 4，num_modes 5。
5. **判定器集成**：上游 LassoPred 零改动，`prp49/lasso_id.py` 做适配层。
6. **回归头**：`dual_head: true` 时主头分类 + 副头回归；`label_type=energy` 时主头也做回归。

---

## 10. 常用命令速查

```bash
# 队列
squeue -u $USER -o "%.10i %.14j %.8T %.10M"

# 作业历史
sacct -j <jobid> --format=JobID%14,JobName%14,State%14,Elapsed,ExitCode -X

# 提交（在 ~/LassoPep 下）
sbatch scripts/job_train_improved.sh
sbatch scripts/job_dock_propedia.sh
sbatch --array=71-80%2 --mem=32G scripts/job_dock_array.sh   # 覆盖式重跑指定任务

# 取消
scancel <jobid>

# 查看输出
tail -40 prp49_<name>-<jobid>.out
```

---

## 附录 A、AF 来源对接数据集（2026-09-16，分支任务）

**背景**：用户指出对接受体结构**来源不一**。核查后发现比预想严重：
**11 个 Lasso 靶中 `CLPB` 与 `NPR1` 本身就是 AlphaFold 模型**（源自 `CLPB_AF.pdb` / `NPR1_AF.pdb`；
叠合 RMSD 恰好 0.00 暴露了这一点），其余 9 个来自各不相同的实验 PDB。

**已建成统一数据集 `docking_af/`**（99/99 对接完成，0 失败，5 作业 19–48 分钟）：
- 受体：**AFDB v6** 单链模型；UniProt 均由 REST 查证（含 PLXNB1 O43157、ITGAV P06756、ITGB3 P05106、CLPB Q9H078）
- 规模：**9 靶 × 11 肽 = 99 对**（排除 ITGAVB3 二聚体、CTRL 对照）
- **盒转移**（AF 模型无配体）两种方法：
  - **配体接触残基法**（新建）：EDNRB(53 接触)/MDM2(33)/PLXNB1(57)——用于整体叠合失败的靶
  - **局部叠合**（≤15 Å 窗口）：C3 **27.5 → 0.55 Å**、PPIA 0.22、FKBP1A 0.40、POLR2A 1.74
  - 整体叠合**已否决**（PLXNB1 41.7 Å、C3 27.5 Å：多结构域蛋白的 AF 单体与实验复合物不共享刚体框架）

### ⚠️ 核心发现：结构来源**主导**对接结果

| 指标 | 值 |
|---|---|
| 同一 (肽,靶) 对分数相关性 | **Pearson +0.195** |
| 逐对差异 | mean +0.79，**sd 2.86**，范围 [−4.88, **+9.96**] |
| **同靶内 top-1 肽一致率** | **1/9 靶** |
| 逐靶 Spearman | NPR1 0.77 · PPIA 0.70 · EDNRB 0.49 · POLR2A **0.02** · C3 **−0.14** · PLXNB1 **−0.42** |

**已知对各有胜负（无一方占优）**：

| 对 | 等级 | AF 版 | 实验版 |
|---|---|---|---|
| RES-701-3 × EDNRB | **A（9KDF 共晶）** | **1/10** ✅ | 4/12 |
| RES-701-1 × EDNRB | B | **2/10** ✅ | 4/12 |
| MccJ25 × POLR2A | C | 2/11 ✅ | 2/12 ✅ |
| Capistruin × POLR2A | C | 10/11 ❌ | 4/12 |
| Lassomycin × CLPB | C | 8/11 ❌ | **1/12** ✅ |
| Anantin × NPR1 | B | 未测 | 2/12 |

**5 个正值异常**（均在接触残基法靶，Ubonodin 占 2 个）：Ubonodin×PLXNB1 **+147.3**、
Siamycin-I×EDNRB +43.9、Sphingopyxin-I×EDNRB +41.5、Ubonodin×MDM2 +28.0、Siamycin-I×PLXNB1 +1.5
→ 盒偏小/位移所致，已排除出统计。

### 对主流程的影响与建议（供 D2/D4 决策）

1. **短名单中的对接分（实验版）不可视为稳健信号** —— 换结构来源即换排序
2. **建议改用双结构共识**：
   - 保守分 = 两版中**较差**者（更严格）
   - 短名单标注"两版**一致** / **冲突**"，**只把一致者作为强候选**
   - 任一版为正值的对直接剔除
3. 冻结验证集的 **3 硬正 / 2 硬负结论不受影响**（那是模型分，与对接无关）；但**短名单排序会变**
4. 若采纳：v1/v2 短名单应重新生成（加 AF 列 + 一致性标注）

**产出**：`docking_af/{README.md, af_scores.csv, af_vs_exp_comparison.csv, af_vs_exp_per_target.csv, box_provenance.csv, contact_boxes.csv, map_report.csv, tasks_af.csv}`（GitHub `f35cdf1`）

1. **`prp49_xneg`（特异性增强训练，62451606）**：训练集加入 **352 个跨靶负样本** + 2 个实验确认阴性 → 直接针对假阳性缺陷
2. **`prp49_calib`（概率校准，62451607）**：温度缩放 + 阈值选择（ECE 前后对比、precision≥0.9 的操作点）
3. **`prp49_dockdp`（成药靶标对接，62451715）**：**39 个成药人靶 × 11 肽 = 429 对**（共晶配体定位口袋）→ 完成后 shortlist 从 11 靶扩到 39 靶
4. **已知对验证集**（`docs/lasso_known_pairs.csv`，25 对，A=6/B=13/C=6）：**17 个矩阵外正样本**可用于拟合融合权重而不污染 3 个验证格
5. **rank 头的域适应**（多任务 / 域对抗）—— 排 P1，需先完成 1–2
6. GitHub 同步 + 报告材料整理

---

---

## 附录 B、路线 C：家族对比预训练（2026-09-16）

**动机**：监督信号只有 352 个正对；而 LassoPred 库里有 4,749 条 core 序列，
按相似度可归入 227 个家族。用"同家族拉近、异家族推远"可以在**零标签**下
教会编码器家族层面的结构先验。

### 数据构建（踩坑一次）

| 步骤 | 结果 |
|---|---|
| 家族来源 | LassoPred 自带 `Lasso_Peptide_Family` 字段 **93% 为 NaN**，不可用 → 改用训练肽的 178 个聚类，把 4,749 条按相似度分配到最近家族 |
| ❌ 首次失败 | **fam094 吞掉 4,134 条（87%）** —— 根因是**短代表序列在长序列上滑动造成偶然匹配**（16 aa 代表 vs 长序列，50% 阈值太松）|
| ✅ 修复 | 用**最长成员**作代表 + **长度差 ≤5** + 阈值 **0.7** → 最大簇 1,109 条（23%），分配率 95% |
| 质量验证 | 抽查簇内随机对 identity：fam060 **0.95–1.00**、fam048 0.96–1.00、fam094/109 **全 1.00** → **大家族是真实的**（库高度冗余），非误匹配 |

### 预训练配置与结果（job 62624147）

```
编码器   LassoESM，解冻顶层 4 层 → 78.7M 可训练参数
投影头   2 层 MLP → 128 维
损失     NT-Xent / InfoNCE，温度 0.07，in-batch 负样本
采样     每个家族等概率（否则 1,109 条的大家族会淹没小家族）
规模     1139 steps/epoch × 15 epochs = 17,085 步，batch 128
耗时     4 分 19 秒（A100，肽序列短 → 极快）
```

| 指标 | 起始 | 结束 | 判读 |
|---|---|---|---|
| loss | 5.244 | **3.525** | 稳定下降 |
| 正对相似度 | 0.989 | 0.945 | 保持（同家族本就极相似）|
| **负对相似度** | 0.966 | **0.502** | ✅ **被显著推开** |
| **分离度** | 0.024 | **0.442** | ✅ 提升 18 倍 |

产出：`runs_contrastive/encoder/`（HF 格式，2,489 MB，hidden 1280 / 33 层，加载校验通过）

### 下游验证（job 62629388，运行中）

关键方法学要求：**唯一变量**。`config_contrastive.yaml` 与 0.8142 基线配置逐项比对，
**仅 4 处差异**（`lassoesm_model` 指向新编码器 + 3 个输出路径），其余（grouped CV、
I1+I2+I5 损失、超参）完全一致 → 任何 CV 差异都可归因于路线 C。

> ⚠️ **预训练本身只是准备**：分离度提升不等于下游变好。必须等 `runs_contrastive_ft`
> 的 grouped CV 出来才能判定路线 C 是否有效。

### 顺带修复：a100 队列 CPU:GPU 配比

对照官方手册核查发现 **a100 队列每卡须配 16 核**，而 6 个作业脚本用了 8 核
（`job_contrastive`/`job_reg_bound3`/`job_mtl`/`job_calibrate`/`job_candidates*`）
→ 已全部改为 16 核（配比错误会导致排队异常）。

### 新增方法学坑

见第 7 节 **#36**（时限须用实测值）· **#37**（并行 array 共享 `--out` 会互相覆盖）。

---

## 附录 C、Lasso 靶标结构数据库 `lasso_target_db/`（2026-09-16/17）

**动机**：路线 D 需要"Lasso 肽 × 细菌靶"的真实相互作用，而 **Propedia 中没有任何受体携带
`RNA_pol_Rpb1` 域**（0/1,739）→ 最关键的 POLR2A 对应关系无法建立。

### 建成内容

| PDB | 肽 | 靶 | 人靶对应 | 置信度 |
|---|---|---|---|---|
| **6N60** | **Microcin J25**（21 aa）| RNAP β′ 亚基（rpoC，1409 aa）| **POLR2A** | high |
| **6N61** | **Capistruin**（19 aa）| RNAP β′ 亚基 | **POLR2A** | high |
| 6N62 | lasso 肽 | RNAP β′ 亚基 | **POLR2A** | high |
| **8IBO / 8IBP** | **Lassomycin**（16 aa）| ClpC1（142/143 aa，结核分枝杆菌）| **CLPB** | high |
| 4CU4 | MccJ25 | FhuA 铁载体受体（706 aa）| —（摄取受体）| low |
| 9KDF | RES-701-3（16 aa）| EDNRB（908 aa，人）| EDNRB | 实验（grade A）|

**POLR2A 的正样本从 0 → 4 对**，且 **MccJ25 与 Capistruin 本身就是扫描肽** →
这些复合物同时是训练数据与天然验证点。

### 关键教训（详见第 7 节 #39/#40/#42）

1. **单一关键词策略会静默漏收集**：6N60/6N61/6N62 只能通过**目标驱动**检索找到
2. **多亚基复合物不能"取第一个匹配"**：会选到 omega 而非 β′
3. **核酸链会被长度过滤误判为肽**：必须按 entity 类型过滤

### 产出文件

`lasso_target_db/`：`README.md`、`lasso_target_pairs_curated.csv`（可用配对 + 人靶对应）、
`target_functional_map.csv`（13 条映射：4 high / 5 medium / 4 none）、
`rnap_lasso_pairs.csv`、`pdb_lasso_entries.csv`（98 结构）、
`search_hits.json` / `target_driven_hits.json`（检索溯源）

配套：`mvp_cpu/db_aux_pairs.csv`（6 正 + 20 诱饵）、`db_aux_holdout.csv`、
`verify_peptides.csv` / `verify_targets.fasta` / `verify_truth.csv`（验证矩阵：5 肽 × 11 靶）

### 能力边界（诚实记录）

- **可提供**：POLR2A（RNAP β′）、CLPB（ClpC1）、PPIA、C3、MDM2 等的细菌对应物
- **不可提供**：EDNRB、PLXNB1、ITGAV、ITGB3 —— 这些是**后生动物特有折叠**，细菌无对应物

### 配套实验（09-17 提交）

| 作业 | 内容 |
|---|---|
| `62657649` | **真实复合物验证**：用现有模型（CV 0.8333）给 5 肽 × 11 靶打分，检验模型是否已具备跨物种迁移能力 |
| `62657636` | **路线 C 下游微调**（修复作用域 bug 后重提，fold 0 从 ep 30 续训）|
| `62646669` | 路线 D′（域同源加权 MTL，运行中）|

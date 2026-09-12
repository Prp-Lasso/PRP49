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

| JobID | 名称 | 内容 | 状态（09-12 09:00 更新）|
|---|---|---|---|
| **62394452** | `prp49_dockp` | 对接批 1 重提（171 个失败任务：18 任务 OOM + 2 超时）| 已改 MAXPAR=4 / mem=190G 重跑 |
| **62363938** | `prp49_dockp2` | 对接批 2（已 1843/1999）| 剩 3 个任务在跑 |
| **62358463** | `prp49_grp` | **无泄漏验证**（按肽分组 CV）| 运行中，**fold0 AUC 0.812 / fold1 AUC 0.849**（还剩 3 折）|
| **62359642** | `prp49_s2` | 跨拓扑迁移（Propedia 8346 对）| fold1 AUC 0.755 / AP 0.761，fold2 进行中 |
| **62394451** | `prp49_regc` | **A+C 亲和力分类**（靶内 z-score + 分箱，2462 对）| 重提（首次因 BatchNorm batch=1 崩溃，已修 drop_last）|
| **62394923** | `prp49_rank` | **成对排序（成药筛选导向）**：16,819 对同靶内偏好样本，5 折 grouped CV | 已提交排队；指标 = pair-acc / 逐靶 Spearman / NDCG@5 / EF@10% |
| 已完成 | `prp49_lassoid2` | Lasso 判定器重训（匹配负样本）| Val F1 **0.7863**；`predict.py` 的 FASTA 注释行 bug 已修并上传 |
| 已停止 | `prp49_reg` | 原亲和力回归（RMSE 6.86，判定不可行）| 已用 A+C 方案取代 |

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
| **grp（无泄漏协议，按肽分组）** | **fold0 0.812 / fold1 0.849** | ✅ **非泄漏确认**：grouped CV 下同样高（还剩 3 折）|
| s2（Propedia 跨拓扑迁移）| fold1 0.755 / AP 0.761 | 迁移学习有信号 |
| regc（A+C 亲和力分类，靶标泛化）| fold0 0.597（仅 1 折，已重提）| 更难的"新靶"任务 |
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

## 5. 下一步待办（按优先级，09-12 09:00 更新）

1. **`grp` 终值**（还剩 3 折）→ 已初步确认非泄漏（fold0 0.812 / fold1 0.849），补齐后写入正式结果
2. **`regc`（A+C 亲和力分类）5 折 AUC** → ≥0.65 说明实验亲和力数据可用于"新靶"预测
3. **对接收尾**：b1 重提（171 对）+ b2 剩 3 → 全量 2,499 完成后重新聚合 `all_scores.csv`
4. **S3 对接增强迁移训练**：配置与数据已就绪（`config_s3.yaml` + `propedia_dock_pairs.csv`），等 GPU 空闲提交 → 对比是否优于 S2
5. **融合分析更新**：用 Propedia 2,161 对对接分 + Lasso 矩阵做 re-ranking 复评
6. **判定器对照集复验**：用随机/打乱/非 Lasso 天然肽验证 F1 0.786（`predict.py` 注释行 bug 已修）
7. **GitHub 同步**（新增：A+C 与 S3 配置、对接聚合表、架构图）
8. **成药面板矩阵**（41 靶 × 13 肽，需先定"成药口袋"盒；AF 模型已下载）
9. **报告/PPT**：验收线 CV ≥ 0.75 ✅ 已过；已知对 top20% ≥ 70% ❌ 未过（当前 50%）

---

## 6. 必须记住的坑（已踩过，勿重蹈）

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

---

## 7. 对接吞吐与资源记录（两轮）

**第一轮（盒 44 Å + exh 4 + 2 核）**：原参数每对需 2+ 小时（2.2 小时仅 398/2499）→ 提速 12–20 倍。
**第二轮（09-12，针对 OOM）**：b1 仍有 18/25 任务 `OUT_OF_MEMORY`（2 个 8 小时超时）→ 再降为 `MAXPAR=4` + `--mem=190G`，失败任务用 `sbatch --array=<job列表>%6` 重提（脚本按 `out/<id>.score` 存在性跳过已完成项）。

---

## 8. 关键设计决策（勿轻易推翻）

1. **评估协议**：主结果用 **grouped CV（按肽序列分组）**；stratified 数字仅作对照（有泄漏）。
2. **硬负样本**：同靶 × 其他家族肽 —— 用于逼模型学"肽-靶兼容性"而非"肽身份"，但因此必须用 grouped CV 评估。
3. **L_align 结构对齐**：用复合物接触图监督 BAN 注意力（`align_loss_weight=0.1`，仅注释对生效）。
4. **对接协议**：AutoDock Vina 1.2.3，刚性配体（大环全冻结，`rigidify.py`），盒=肽自身坐标 + pad，exhaustiveness 4，num_modes 5。
5. **判定器集成**：上游 LassoPred 零改动，`prp49/lasso_id.py` 做适配层。
6. **回归头**：`dual_head: true` 时主头分类 + 副头回归；`label_type=energy` 时主头也做回归。

---

## 9. 常用命令速查

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

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

| JobID | 名称 | 内容 | 预期 |
|---|---|---|---|
| **62363937** | `prp49_dockp` | 对接批 1 剩余（500 对中已完成 56）| 1–1.5 小时 |
| **62363938** | `prp49_dockp2` | 对接批 2 剩余（1999 对中已完成 342）| 3–5 小时 |
| **62358463** | `prp49_grp` | **无泄漏验证**（按肽分组 CV，3248 行）| 5–7 小时 |
| **62359642** | `prp49_s2` | 跨拓扑迁移（Propedia 8346 对）| 数小时 |
| **62363934** | `prp49_reg` | **亲和力回归**（4321 对实验标签，已修 BCE bug）| 2–4 小时 |
| 待确认 | `prp49_lassoid2` | Lasso 判定器重训（匹配负样本版）| **训练已完成**：Val F1 **0.7863**；best_model.pt 已落盘（09-11 19:17）；仅后续 predict 阶段因 FASTA 注释行报错（可选修）|

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
| **I1+I2+I5 全套改进** | **0.8145 ± 0.017** | ⚠️ **疑似肽级泄漏，待 grp 验证** |
| grp 首轮（无泄漏协议，训练中）| 0.761 | 中途数字，待终值 |
| Lasso 判定器（原版，长度捷径）| F1 **0.973** / AUC 0.998 | ❌ 负样本长度分布不重叠 → 虚高 |
| **Lasso 判定器（匹配负样本，诚实版）** | **Val F1 0.7863** | ✅ best_model.pt 已是此版本（CD-HIT 去重：正 1697 / 负 4414；训练 4888 对）|

**对接矩阵**：11 肽 × 10 靶 = 110 对完成（`docking/dock_matrix_full.csv`）
- 阳性对照 CTRL_9KDF redocking = **−12.2**（盒放大后合理；放大前是 +104.8 假信号）
- 模型分 vs 对接分 Spearman **ρ = −0.091**（不相关 → 信息互补）
- 融合 `z(model)+w·z(−dock)`：已知对平均排名 **4.43 → 3.57**（w=1.0）；前 20% 命中 4/8 = 50%（**未达验收 70%**）
- 主失败案例：**PB1m7×PLXNB1 排 11/11**（人工嫁接肽，模型盲点）、MccJ25×POLR2A 排 10/11

---

## 4. 数据资产（服务器 `~/LassoPep/`）

| 路径 | 内容 |
|---|---|
| `mvp_cpu/train_pairs_hard.csv` | 训练集（正 312 + 硬负样本），improved 用 |
| `mvp_cpu/affinity_pairs.csv` | **4321 对实验亲和力**（Ki 1615 / Kd 545 / IC50 2161；564 靶；标签列 `energy` = pAffinity 3–12）|
| `mvp_cpu/scan_peptides_ext.csv` | 13 条扫描肽 |
| `mvp_cpu/scan_targets.fasta` | 11 条靶链 |
| `data/lassopred_database.csv` | LassoPred 库（拓扑注释来源，4749 条；**注意路径是 `data/` 不是根目录**）|
| `lasso_id_data/` | 判定器数据（正 4749 / 负 20000 / 匹配负 10559）|
| `LassoPeptideClassifier/` | 判定器代码 + `checkpoints/`（待 `lassoid2` 落盘）|
| `docking/` | Lasso 矩阵（110 对）+ `tasks.csv` + `pep/`、`rec/`、`pdbqt/`、`out/` |
| `docking_propedia/` | 批 1（500 对实验复合物）|
| `docking_propedia_b2/` | 批 2（1999 对）|
| `docking_drugpanel/af/` | 41 个成药人靶 AF 模型（备用）|
| `PRP49/runs*/checkpoints/fold*_best.pt` | 各训练权重 |
| `results/scan_matrix_ext.csv` | 模型打分矩阵（13 肽）|

**本机对应路径**：`D:\deepseek_harness\prp49\`（`mvp_cpu/`、`docking/`、`results/`、`scratch/`（分析脚本）、`affinity_data/`）

---

## 5. 下一步待办（按优先级）

1. **对接完成**（2499 对）→ 聚合出矩阵 → 转回归标签（`energy`）→ 提交**回归训练 v2**
2. **`grp` 终值** → 判定 0.8145 是否为泄漏；若跌回 ~0.6，需重设计硬负样本
3. **判定器对照集复验**：用随机序列 / 打乱序列 / 非 Lasso 天然肽测判定器，确认不是长度捷径
4. **`reg` 结果**（Spearman / RMSE）→ 决定是否扩大亲和力数据
5. **融合分析更新**（用 2499 对对接分重做 re-ranking）
6. **GitHub 同步**（`github_repo/` → https://github.com/Prp-Lasso/PRP49，SSH 密钥已配好，`git push` 前先 `git add`）
7. **成药面板矩阵**（41 靶 × 13 肽，需先定"成药口袋"盒；AF 模型已下载）
8. **报告/PPT**：把当前成果整理成汇报（验收线：CV ≥ 0.75 ✓ 已过；已知对 top20% ≥ 70% ✗ 未过）

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

---

## 7. 对接吞吐优化记录（本轮）

原参数每对需 2+ 小时（2.2 小时仅 398/2499）。优化：
- 盒边长上限 **60 → 44 Å**（b1 均体积 216k → 52.8k Å³；b2 → 37.8k Å³）
- `--exhaustiveness 8 → 4`
- 每任务 `--cpu 1 → 2`，并发 b1 `%4→%6`、b2 `%8→%12`

→ 综合提速 12–20 倍。已完成的 398 个结果保留（脚本按 `out/<id>.score` 存在性跳过）。

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

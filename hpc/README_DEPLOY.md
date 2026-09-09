# PRP49 交我算部署手册（账号到手后按序执行）

> 依据官方手册 https://docs.hpc.sjtu.edu.cn 整理。所有占位符 `YOU` 需替换为真实 jAccount 账号、邮箱。
> 目标集群：**思源一号**（a100 队列，A100 40GB）；π 2.0 仅作对照/备用。

## 黄金规则（违反会封号）

1. **登录节点严禁跑计算**——一切训练/评估必须 `sbatch` 到计算节点；登录节点只做编辑、下载、提交。
2. 计算节点无公网直连，下载必须走代理（思源一号 `proxy2.pi.sjtu.edu.cn:3128`；π 2.0 `proxy.pi.sjtu.edu.cn:3004`）。本项目模型已本地打包，**无需下载**。
3. 文件/目录名只用英文+数字。
4. 无 sudo；软件用 module + conda。
5. `$SCRATCH` 每 3 个月清理，重要产物及时 `cp` 回 `$HOME`（`/lustre`）。
6. a100 配比：**每卡 16 CPU 核**（`--cpus-per-task=16`），写错会长期 PD 或浪费机时。

---

## Step 0 — 登录与集群体检（登录节点）

```bash
# 思源一号（主集群）
ssh YOU@sylogin1.hpc.sjtu.edu.cn
# π 2.0（备用/对照）
ssh YOU@pilogin.hpc.sjtu.edu.cn
```

登录后体检（只读命令，不算计算）：

```bash
sinfo --partition=a100            # a100 节点状态，idle=可用
scontrol show partition a100      # 确认 MaxTime（默认 3 天；debuga100 时限以这里为准）
squeue -u $USER                   # 自己无遗留作业
echo $SCRATCH                     # 确认 SCRATCH 路径（/scratch/home/...）
module avail miniconda            # 确认 miniconda3/4.10.3 存在
```

## Step 1 — 上传数据（本机 → 思源一号）

在本机（`hpc_deploy/upload/`）：

```powershell
# 1) 打包（生成 prp49_bundle.tar.gz + SHA256）
.\make_bundle.ps1
# 2) 上传（约 2.7 GB；scp 走校园网约 10-30 分钟）
scp prp49_bundle.tar.gz YOU@sylogin1.hpc.sjtu.edu.cn:~/
```

在集群登录节点：

```bash
sha256sum ~/prp49_bundle.tar.gz    # 与本机 SHA256 比对
mkdir -p ~/prp49
tar -xzf ~/prp49_bundle.tar.gz -C ~/prp49
ls ~/prp49                          # 应见：PRP49/ LassoESM_hf/ esm2_35M_hf/ mvp_cpu/ data/ scripts/ hpc_deploy/ docs/ TUnA/
```

> 大文件断点续传备选：`rsync -avhP --partial prp49_bundle.tar.gz YOU@sylogin1.hpc.sjtu.edu.cn:~/`
> 若 tar.gz 超过网络容忍，可只传模型目录：`rsync -avhP --partial LassoESM_hf esm2_35M_hf YOU@sylogin1.hpc.sjtu.edu.cn:~/prp49/`

## Step 2 — conda 环境（登录节点只建环境；装包可放交互节点）

```bash
module load miniconda3/4.10.3
conda create -n prp49 python=3.10 -y
source activate prp49
cd ~/prp49
conda env update -n prp49 -f hpc_deploy/environment.yml
```

> environment.yml / requirements.txt 在 bundle 的 `hpc_deploy/` 下（bundler 会带）。若没带，用下面手动装：

```bash
source activate prp49
pip install numpy pandas scikit-learn pyyaml tqdm
pip install transformers==4.40.2 huggingface_hub
# PyTorch CUDA 版（A100 用 cu121；走代理或清华镜像二选一）
export https_proxy=http://proxy2.pi.sjtu.edu.cn:3128 http_proxy=http://proxy2.pi.sjtu.edu.cn:3128
pip install torch==2.3.1 --index-url https://download.pytorch.org/whl/cu121
# 若代理下载慢：pip install torch==2.3.1 -i https://pypi.tuna.tsinghua.edu.cn/simple（注意 pypi 默认版已带 CUDA 运行时）
```

验证环境（在交互节点，见下）：
```bash
python -c "import torch; print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

## Step 3 — 冒烟测试（debuga100，≤20 分钟）

```bash
cd ~/prp49/scripts
sbatch job_smoke_debuga100.sh
squeue -u $USER     # 等到 R → 结束
cat prp49_smoke-*.out | tail -30   # 应见 torch cuda: True + SMOKE_DONE
```

目的：验证 GPU、环境、数据路径（config_smoke.yaml 用 8M 模型，5G 虚拟卡足够）。

## Step 4 — 正式训练（a100 单卡起步）

```bash
cd ~/prp49/scripts
# 改好作业脚本里的 --mail-user=YOU@sjtu.edu.cn（所有 job_*.sh）
sbatch job_train_a100.sh
squeue -u $USER          # PD=排队正常；超过 10 分钟仍 PD → 检查配比/配额
```

**首轮执行范围**（与 docs/超算首轮执行与文档差距清单.md 一致）：

- **S1 预热（Bernett L0，~1–2 GPU·h）**：`python -m prp49.train --config config_warmup.yaml --device cuda`（TUnA Bernett 5 万对已随 bundle 在 `mvp_cpu/bernett_pairs.csv`，双 ESM-2 35M 冻结只训 BAN）→ 产出 warmup 检查点作 S2/S3 初始化；
- **L2 微调基线**：`job_train_a100.sh` 用 config.yaml（expanded 312 正），与 8M MVP（分层 0.663/分组 0.641）同协议对标；
- 含增强集的变体训练（可选）：`sed 's/expanded_pairs.csv/train_pairs.csv/' config.yaml > config_aug.yaml` 后把作业脚本里 config 换成 config_aug.yaml 提交（train_pairs = 1248 对 + 40 保守突变变体）；
- Propedia（L1 跨拓扑环肽-人靶）与 TUnA 对照复现为第二轮（见差距清单 A1/A4）。

- 单次 5 折 CV 估算 100–180 GPU·h（3090 口径全流程值）；A100 40GB 约 1.5–2 倍速，单次 CV 预计 1–2 天。`--time=2-00:00:00` 不够就提前 1 个工作日邮件申请延长（≤14 天）。
- 多 seed / 多配置并行：`sbatch job_train_array.sh`（`--array=1-4%4`，每个 task 独立 config、独立输出目录）。
- checkpoint 在 `~/prp49/PRP49/runs/checkpoints/fold{0..4}_best.pt`（软链到 `$SCRATCH/prp49_runs`，快但 3 个月清理）。
- 训练完 `cp $SCRATCH/prp49_runs/cv_summary.json ~/prp49/results/` 留档。

## Step 5 — 评估 + 消融 + 回溯 + 10×10 扫描

```bash
cd ~/prp49/scripts
sbatch job_eval.sh runs/checkpoints/fold0_best.pt     # held-out AUC + 环消融(Wilcoxon/80%判据)
sbatch job_scan.sh runs/checkpoints/fold0_best.pt     # 输出 ../results/scan_matrix.csv
# 回溯验证（已知对排名百分位，验收 ≥70% 进靶标前 20%）：训练完成后回本地或交互节点跑
python -m prp49.retrospective --config config.yaml --ckpt runs/checkpoints/fold0_best.pt \
    --retro ../mvp_cpu/retro_pairs.csv --peptides ../mvp_cpu/scan_peptides.csv --out retro_rank.csv
```

> `retro_pairs.csv`（列 pep_id,pep_seq,prot_seq）为留出的已知互作对清单（9KDF/7VF3/6N60/6N61/8IBO 等，训练时未含的 fold）；扫描肽池作排名池，可加 `--neg_csv` 扩充同靶负样本池。
> 全蛋白组粗筛：`../data/human_proteome.fasta`（20,416 条）可直接作为 `job_scan.sh` 的 `--targets`（改脚本或另提交），输出 10 肽 × 2 万蛋白矩阵约几十 MB。

## Step 6 — 结果回传与清理

```bash
# 集群 → 本机
scp YOU@sylogin1.hpc.sjtu.edu.cn:~/prp49/results/scan_matrix.csv .
scp -r YOU@sylogin1.hpc.sjtu.edu.cn:~/prp49/PRP49/runs .
# SCRATCH 空间管理：训练中间产物 3 个月自动清理；要保留的先 cp 回 /lustre
```

## 故障速查

| 现象 | 处理 |
|---|---|
| 作业长期 PD | a100 配比必须 `--cpus-per-task=16 --gres=gpu:1`；队列忙则等待 |
| `node_fail` | 节点故障，直接 `sbatch` 重投，机时自动返还 |
| 计算节点下载失败 | 没配代理（本项目模型已本地打包，无需下载） |
| 登录节点跑程序被 kill | 违反黄金规则 1，等解封（30–120 min）后 sbatch |
| `ModuleNotFoundError: prp49` | 没在 `~/prp49/PRP49` 目录下运行（作业脚本已 cd，手动跑注意 cwd） |
| OOM | LassoESM 650M + ESM-2 35M 单卡 batch 16 理论 OK；报 OOM 就 batch 降到 8 或 `lassoesm_unfreeze: 0` |
| 需要延长时限 | 提前 1 个工作日邮件 hpc 邮箱：用户名 + jobID |

## 支持

- 文档 https://docs.hpc.sjtu.edu.cn ｜ AI 助手 https://chat.hpc.sjtu.edu.cn ｜ 门户 https://hpc.sjtu.edu.cn
- 论文致谢（英文）：The computations in this paper were run on the Siyuan-1 cluster supported by the Center for High Performance Computing at Shanghai Jiao Tong University.

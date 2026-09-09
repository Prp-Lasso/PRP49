#!/bin/bash
# ============================================================
# PRP49 正式训练：思源一号 a100 队列，单卡
# 用法：cd ~/prp49/scripts && sbatch job_train_a100.sh
# 配比：a100 队列每卡配 16 CPU 核（官方规定，勿改）
# ============================================================
#SBATCH --job-name=prp49_train
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00          # 单次 CV 预估 100-180 GPU·h，超时可提前 1 个工作日邮件申请延长（≤14 天）
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOU@sjtu.edu.cn   # ← 改成你的邮箱

set -e

# 模型/数据全部本地化，禁止联网（避免卡在 HF 探测超时）
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

PRP49_ROOT=$HOME/prp49
cd $PRP49_ROOT/PRP49

# checkpoint / 日志写 $SCRATCH（全闪存，快；3 个月清理），软链回代码目录
mkdir -p $SCRATCH/prp49_runs
ln -sfn $SCRATCH/prp49_runs ./runs

# 环境（思源一号 miniconda 版本）
module load miniconda3/4.10.3
source activate prp49

echo "=== node: $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv
echo "=== torch: $(python -c 'import torch; print(torch.__version__, torch.cuda.is_available())') ==="

python -m prp49.train --config config.yaml --device cuda

echo "TRAIN_DONE"

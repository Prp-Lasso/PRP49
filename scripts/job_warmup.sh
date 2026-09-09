#!/bin/bash
# ============================================================
# PRP49 S1 预热：Bernett 通用 PPI（5 万对），a100 单卡 ~1-2 小时
# 双 ESM-2 35M 全冻结，只训练 BAN + 分类头
# 产出：runs_warmup/checkpoints/（作 S2/S3 初始化，后续可接课程学习）
# ============================================================
#SBATCH --job-name=prp49_warmup
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=00-03:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOU@sjtu.edu.cn   # ← 改成你的邮箱

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

PRP49_ROOT=$HOME/prp49
cd $PRP49_ROOT/PRP49

mkdir -p $SCRATCH/prp49_runs
ln -sfn $SCRATCH/prp49_runs ./runs_warmup

module load miniconda3/4.10.3
source activate prp49

echo "=== node: $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv
echo "=== pairs: $(wc -l ../mvp_cpu/bernett_pairs.csv) ==="

python -m prp49.train --config config_warmup.yaml --device cuda

echo "WARMUP_DONE"

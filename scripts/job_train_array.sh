#!/bin/bash
# ============================================================
# PRP49 多 seed 并行：--array=1-4%4 → 4 个独立单卡作业（每作业 1 卡 16 核）
# 每个 array task 生成独立 config（不同 seed / 输出目录），互不干扰
# %4 = 最多同时跑 4 个（按账号 GPU 配额调整，配额不足会自动排队）
# ============================================================
#SBATCH --job-name=prp49_seed
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --array=1-4%4
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOU@sjtu.edu.cn   # ← 改成你的邮箱

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

PRP49_ROOT=$HOME/prp49
cd $PRP49_ROOT/PRP49

mkdir -p $SCRATCH/prp49_runs
ln -sfn $SCRATCH/prp49_runs ./runs

module load miniconda3/4.10.3
source activate prp49

SEED=$((40 + SLURM_ARRAY_TASK_ID))
CONF=config_seed${SEED}.yaml
sed -e "s/^  seed: .*/  seed: ${SEED}/" \
    -e "s|^  output_dir: .*|  output_dir: ./runs_seed${SEED}|" \
    -e "s|^  checkpoint_dir: .*|  checkpoint_dir: ./runs_seed${SEED}/checkpoints|" \
    config.yaml > ${CONF}

echo "=== task ${SLURM_ARRAY_TASK_ID} seed=${SEED} conf=${CONF} ==="
python -m prp49.train --config ${CONF} --device cuda
echo "SEED_DONE ${SEED}"

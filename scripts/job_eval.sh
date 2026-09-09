#!/bin/bash
# ============================================================
# PRP49 评估 + 环依赖消融：a100 单卡
# 用法：sbatch job_eval.sh <checkpoint.pt>
#   ckpt 来自训练作业：runs/checkpoints/fold{0..4}_best.pt（软链到 $SCRATCH）
# ============================================================
#SBATCH --job-name=prp49_eval
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=00-02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

CKPT=${1:?usage: sbatch job_eval.sh <runs/checkpoints/fold0_best.pt>}

cd $HOME/prp49/PRP49
module load miniconda3/4.10.3
source activate prp49

python -m prp49.eval --config config.yaml --ckpt ${CKPT} --device cuda
echo "EVAL_DONE"

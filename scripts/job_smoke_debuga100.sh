#!/bin/bash
# ============================================================
# PRP49 冒烟测试：思源一号 debuga100 调试队列
# 目的：验证 GPU、conda 环境、数据路径、代码入口全部可用
# 注意：debuga100 是虚拟小卡（每卡 5G 显存），只能用 config_smoke.yaml
#      （双 8M ESM 冻结编码器，显存 <1G）
# 每人同时最多 1 个 debug 作业；时限以 scontrol show partition debuga100 为准
# ============================================================
#SBATCH --job-name=prp49_smoke
#SBATCH --partition=debuga100
#SBATCH --qos=debug                  # debuga100 必加
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=00:20:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

source /usr/share/lmod/lmod/init/profile
module load miniconda3/4.10.3
source activate prp49

PRP49_ROOT=$HOME/LassoPep
cd $PRP49_ROOT/PRP49

echo "=== host: $(hostname) ==="
nvidia-smi || true
echo "=== python: $(python --version 2>&1) ==="
echo "=== torch cuda: $(python -c 'import torch; print(torch.__version__, torch.cuda.is_available())') ==="
echo "=== transformers: $(python -c 'import transformers; print(transformers.__version__)') ==="

# 数据路径检查（config_smoke.yaml 依赖这些文件）
ls -la ../LassoESM_hf/model.safetensors ../esm2_8M_hf/model.safetensors ../mvp_cpu/pairs.csv ../mvp_cpu/contacts/contacts.json

python -m prp49.train --config config_smoke.yaml --device cuda

echo "SMOKE_DONE"

#!/bin/bash
# ============================================================
# 备选：计算节点下载模型作业（仅当本地打包上传不可用时用）
# 注意：本作业在【π 2.0 登录节点 pilogin.hpc.sjtu.edu.cn】提交，
#      代理是 π 2.0 的 proxy.pi.sjtu.edu.cn:3004（与思源一号不同！）
# 主方案：本地已打包 esm2_35M_hf + LassoESM_hf 随 bundle 上传，无需此脚本
# ============================================================
#SBATCH --job-name=prp49_dl
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00-30:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -e
export http_proxy=http://proxy.pi.sjtu.edu.cn:3004/
export https_proxy=http://proxy.pi.sjtu.edu.cn:3004/
export no_proxy=puppet

cd $HOME/prp49
module load miniconda3/4.8.2
source activate prp49

python - <<'EOF'
from huggingface_hub import snapshot_download
p = snapshot_download('facebook/esm2_t12_35M_UR50D', local_dir='./esm2_35M_hf')
print('DONE', p)
EOF

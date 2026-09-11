"""Create prp49 conda env on login node: python 3.10 + torch cu121 + transformers."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

INIT = 'source /usr/share/lmod/lmod/init/profile; module load miniconda3/4.10.3; '
TSINGHUA = '-i https://pypi.tuna.tsinghua.edu.cn/simple'

cli = connect()

print('=== [1/4] conda create prp49 (py3.10) ===', flush=True)
st, out, err = run(cli, INIT + 'conda create -n prp49 python=3.10 -y 2>&1 | tail -3', timeout=600)
print(out, err, flush=True)

print('=== [2/4] pip base packages ===', flush=True)
st, out, err = run(cli, INIT + 'source activate prp49 && pip install ' + TSINGHUA +
                   ' numpy pandas scikit-learn pyyaml tqdm transformers==4.40.2 huggingface_hub 2>&1 | tail -3',
                   timeout=900)
print(out, err, flush=True)

print('=== [3/4] pip torch (cu121 default wheel, ~2.4GB) ===', flush=True)
st, out, err = run(cli, INIT + 'source activate prp49 && pip install ' + TSINGHUA +
                   ' torch==2.3.1 2>&1 | tail -3', timeout=2400)
print(out, err, flush=True)

print('=== [4/4] verify ===', flush=True)
st, out, err = run(cli, INIT + 'source activate prp49 && python -c "'
                   'import torch, transformers, sklearn, pandas, numpy; '
                   'print(\'torch\', torch.__version__, \'cuda_build\', torch.version.cuda, '
                   '\'transformers\', transformers.__version__)" 2>&1', timeout=300)
print(out, err, flush=True)

cli.close()
print('SETUP_ENV DONE', flush=True)

"""Upload improved-run assets (code+config+data+job) and submit."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

cli = connect()
st, out, err = run(cli, 'echo $HOME')
home = out.strip()

for local, remote in [
    (r'D:\deepseek_harness\prp49\PRP49\prp49\train.py', 'PRP49/prp49/train.py'),
    (r'D:\deepseek_harness\prp49\PRP49\prp49\losses.py', 'PRP49/prp49/losses.py'),
    (r'D:\deepseek_harness\prp49\PRP49\config_improved.yaml', 'PRP49/config_improved.yaml'),
    (r'D:\deepseek_harness\prp49\mvp_cpu\train_pairs_hard.csv', 'mvp_cpu/train_pairs_hard.csv'),
    (r'D:\deepseek_harness\prp49\hpc_deploy\scripts\job_train_improved.sh', 'scripts/job_train_improved.sh'),
]:
    put_file(cli, local, f'{home}/LassoPep/{remote}', progress=False)

st, out, err = run(cli, 'cd ~/LassoPep && sed -i "s/\\r$//" mvp_cpu/train_pairs_hard.csv; '
                        'wc -l mvp_cpu/train_pairs_hard.csv; chmod +x scripts/job_train_improved.sh; '
                        'sbatch scripts/job_train_improved.sh')
print('improved submit:', out.strip(), err.strip())
st, out, err = run(cli, 'squeue -u $USER -o "%.10i %.14j %.8T %.10M" | head -8')
print(out)
cli.close()

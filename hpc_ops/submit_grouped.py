"""Submit leakage-safe grouped-CV verification run."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

cli = connect()
st, out, err = run(cli, 'echo $HOME')
home = out.strip()
put_file(cli, r'D:\deepseek_harness\prp49\PRP49\config_improved_grouped.yaml',
         f'{home}/LassoPep/PRP49/config_improved_grouped.yaml', progress=False)
put_file(cli, r'D:\deepseek_harness\prp49\hpc_deploy\scripts\job_train_grouped.sh',
         f'{home}/LassoPep/scripts/job_train_grouped.sh', progress=False)
st, out, err = run(cli, 'cd ~/LassoPep && chmod +x scripts/job_train_grouped.sh; '
                        'sbatch scripts/job_train_grouped.sh; sleep 3; '
                        'squeue -u $USER -o "%.10i %.14j %.8T %.10M"')
print(out, err)
cli.close()

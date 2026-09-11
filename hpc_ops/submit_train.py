"""Submit warmup + main training on a100 and report queue status."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

cli = connect()
print('=== submit ===', flush=True)
st, out, err = run(cli, 'cd ~/LassoPep && sbatch scripts/job_warmup.sh; sbatch scripts/job_train_a100.sh')
print(out, err, flush=True)

print('=== queue ===', flush=True)
st, out, err = run(cli, 'squeue -u $USER -o "%.10i %.14j %.8T %.10M %.6D %.20R"')
print(out, flush=True)

print('=== a100 availability ===', flush=True)
st, out, err = run(cli, 'sinfo --partition=a100 -o "%.10P %.6t %.10N %.6D %.20E" | head -8')
print(out, flush=True)
cli.close()
print('SUBMIT_TRAIN DONE', flush=True)

"""Resubmit warmup with 20k subset + shorter epochs (fits 3h slot)."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

cli = connect()
st, out, err = run(cli, 'echo $HOME')
home = out.strip()

# make 20k subset on server (file op, not compute)
st, out, err = run(cli, 'cd ~/LassoPep/mvp_cpu && head -20001 bernett_pairs.csv > bernett_pairs_20k.csv && wc -l bernett_pairs_20k.csv')
print('subset:', out.strip(), err.strip())

put_file(cli, r'D:\deepseek_harness\prp49\PRP49\config_warmup.yaml',
         f'{home}/LassoPep/PRP49/config_warmup.yaml', progress=False)

# cancel old warmup, resubmit
st, out, err = run(cli, 'scancel 62206423; sleep 2; cd ~/LassoPep && sbatch scripts/job_warmup.sh')
print('resubmit:', out.strip(), err.strip())

st, out, err = run(cli, 'squeue -u $USER -o "%.10i %.14j %.8T %.10M %.20R"')
print(out)
cli.close()
print('WARMUP_RESUBMIT DONE')

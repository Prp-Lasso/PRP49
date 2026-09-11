"""Probe: login, verify lassopep dir, cluster health. Read-only."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

cli = connect()
st, out, err = run(cli, 'hostname; echo "---"; whoami; echo "---"; echo $HOME')
print('login:', st, out.strip(), err.strip())

st, out, err = run(cli, 'ls -la ~/lassopep 2>&1 | head -30')
print('--- lassopep dir ---')
print(out)

st, out, err = run(cli, 'sinfo --partition=a100 2>&1 | head -5; echo "==="; squeue -u $USER 2>&1 | head -10; echo "==="; df -h $HOME 2>&1 | tail -2')
print('--- cluster state ---')
print(out)
cli.close()
print('PROBE DONE')

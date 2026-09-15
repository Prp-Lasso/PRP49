"""Sync server-side job scripts back to local, re-run preflight, push to GitHub."""
import os
import subprocess
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

ROOT = r'D:\deepseek_harness\prp49'
cli = connect()
st, h, _ = run(cli, 'echo $HOME')
h = h.strip()

# server is the source of truth for job scripts (edited in place during fixes)
for f in ['job_reg_bound2.sh', 'job_train_prop.sh', 'job_train_family.sh', 'job_mtl.sh']:
    try:
        sftp = cli.open_sftp()
        sftp.get(f'{h}/LassoPep/scripts/{f}', os.path.join(ROOT, 'hpc_deploy', 'scripts', f))
        sftp.close()
        print(f'synced {f} from server')
    except Exception as exc:
        print(f'{f}: {exc}')
cli.close()

# sanitize the local copy of the data file so a future upload cannot reintroduce CRLF
target = os.path.join(ROOT, 'mvp_cpu', 'affinity_pairs.csv')
raw = open(target, 'rb').read()
n = raw.count(b'\r\n')
if n:
    open(target, 'wb').write(raw.replace(b'\r\n', b'\n'))
    print(f'sanitized {target}: removed {n} CRLF')

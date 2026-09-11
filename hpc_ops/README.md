# HPC remote-ops helpers (SJTU Siyuan-1)

Thin paramiko wrappers used to drive the cluster from a workstation:
upload, set up conda envs, submit/monitor Slurm jobs, fetch results.

Golden rules enforced by these scripts (see docs/超算首轮执行与文档差距清单.md):
- heavy work only via `sbatch`; the login node is used for transfers and submissions
- compute nodes have no internet; set http(s)_proxy in job scripts when downloading
- stay strictly inside the user's own directory (BASE)

Files
- `opslib_template.py` - connection helper (credentials from env vars)
- `probe.py`           - login check, partition/queue health
- `setup_env.py`       - build the `prp49` conda env (torch cu121 + transformers)
- `submit_*.py`        - submit training / scan / docking arrays
- `wait_dock3.py`, `poll_bigbox2.py` - poll arrays and aggregate scores
- `agg_full_dock.py`, `get_scan.py`  - collect matrices back to local

Usage: `python probe.py` after exporting SJTU_HPC_USER / SJTU_HPC_PASSWORD.

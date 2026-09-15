"""Static verification of the checkpoint wiring in every training script."""
import ast
import os
import re

ROOT = r'D:\deepseek_harness\prp49\PRP49\prp49'
NAMES = ['ckpt_dir', 'ckpt_every', 'start_ep']

for f in ['train.py', 'train_reg_bound.py', 'train_mtl.py', 'train_pair.py']:
    p = os.path.join(ROOT, f)
    s = open(p, encoding='utf-8').read()
    try:
        ast.parse(s)
        syn = 'syntax OK'
    except SyntaxError as exc:
        print(f'{f}: SYNTAX ERROR {exc}')
        continue
    imp = 'from .ckpt import' in s
    print(f'{f:22s} {syn} | import={imp} | '
          f'maybe_resume={s.count("maybe_resume(")} save_resume={s.count("save_resume(")} '
          f'save_best={s.count("save_best(")}')
    for name in NAMES:
        used = len(re.findall(rf'\b{name}\b', s))
        defined = len(re.findall(rf'\b{name}\s*=', s))
        if used and not defined:
            print(f'    !! {name}: used {used}x but NEVER defined -> NameError at runtime')
        elif used:
            print(f'    {name}: used {used}x, defined {defined}x OK')

# sanity: does the ckpt helper itself behave as intended?
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\PRP49')
os.chdir(r'D:\deepseek_harness\prp49\PRP49')
print('\n--- ckpt helper smoke ---')
import tempfile
import torch
import torch.nn as nn
from prp49.ckpt import ckpt_paths, maybe_resume, save_atomic, save_best, save_resume

tmp = tempfile.mkdtemp()
model = nn.Linear(4, 2)
opt = torch.optim.Adam(model.parameters())
save_best(tmp, 0, model.state_dict(), {'auc': 0.8})
save_resume(tmp, 0, model, opt, epoch=9, best=0.8)
latest, best = ckpt_paths(tmp, 0)
print('files:', sorted(os.listdir(tmp)))
m2, o2 = nn.Linear(4, 2), torch.optim.Adam(model.parameters())
start, b = maybe_resume(tmp, 0, m2, o2, device='cpu', verbose=True)
print(f'resume -> start_ep={start} best={b} (expect start_ep=10 best=0.8)')
stray = [x for x in os.listdir(tmp) if x.endswith('.tmp')]
print('leftover .tmp files:', stray or 'none')

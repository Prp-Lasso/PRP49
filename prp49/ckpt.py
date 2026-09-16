"""Checkpoint utilities: periodic saving + resume, safe against wall-clock kills.

Why this exists
---------------
Every training script used to save only at the END of a fold. Three jobs were
killed by their time limit (regb serial, regb2 8h, MTL 6h) and in each case all
in-memory progress - including the best epoch - was lost. Two rules fix that:

  1. SAVE PERIODICALLY. Every `ckpt_every` epochs the resume checkpoint is
     overwritten, so a killed job always leaves a usable state behind.
  2. WRITE ATOMICALLY. Data goes to a .tmp file that is renamed into place, so a
     kill during the write cannot leave a truncated checkpoint.

Two files per fold, both full-model for compatibility:
  fold{i}_latest.pt - rolling state + optimiser + epoch   (for resume)
  fold{i}_best.pt   - best validation state               (for delivery)
"""
import os

import torch


def ckpt_paths(ckpt_dir, fold):
    return (os.path.join(ckpt_dir, f'fold{fold}_latest.pt'),
            os.path.join(ckpt_dir, f'fold{fold}_best.pt'))


def save_atomic(obj, path):
    """Write via a temp file + rename so an interrupted write cannot corrupt it."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp'
    torch.save(obj, tmp)
    os.replace(tmp, path)


def _trainable_state(model):
    """State dict restricted to parameters that actually train.

    The encoders are largely frozen (LassoESM with only the top layer unfrozen,
    ESM-2 fine-tuned but small), so storing the full 2.75 GB state on every rolling
    save wastes IO - and a rolling file that also carries Adam moments reached
    3.18 GB, slowing training on GPFS. Keeping only trainable tensors shrinks the
    rolling checkpoint by roughly an order of magnitude.
    """
    names = {n for n, p in model.named_parameters() if p.requires_grad}
    return {k: v for k, v in model.state_dict().items() if k in names}


def save_resume(ckpt_dir, fold, model, optimizer, epoch, best, extra=None):
    """Rolling checkpoint: trainable weights + epoch, without optimiser moments.

    Optimiser state is deliberately dropped - it would double the size, and the
    Adam moments are cheap to rebuild (a few epochs of slightly slower convergence
    after a resume, in exchange for a checkpoint small enough to write often).
    """
    latest, _ = ckpt_paths(ckpt_dir, fold)
    save_atomic({'epoch': int(epoch), 'state_dict': _trainable_state(model),
                 'optimizer': None, 'best': float(best), 'extra': extra or {},
                 'partial_state': True}, latest)
    return latest


def save_best(ckpt_dir, fold, state_dict, metrics=None):
    _, best = ckpt_paths(ckpt_dir, fold)
    save_atomic({'state_dict': state_dict, 'metrics': metrics or {}}, best)
    return best


def maybe_resume(ckpt_dir, fold, model, optimizer=None, device='cuda', verbose=True):
    """Load fold{i}_latest.pt if present. Returns (start_epoch, best)."""
    latest, _ = ckpt_paths(ckpt_dir, fold)
    if not os.path.exists(latest):
        return 0, -1.0
    try:
        ck = torch.load(latest, map_location=device)
    except Exception as exc:                       # corrupt file -> start clean
        if verbose:
            print(f'  [resume] fold {fold}: cannot read {latest} ({exc}); starting fresh',
                  flush=True)
        return 0, -1.0
    try:
        sd = ck['state_dict']
        if ck.get('partial_state'):
            # rolling checkpoints now hold only trainable tensors: merge them into
            # the freshly built model instead of replacing the whole state dict
            full = model.state_dict()
            full.update({k: v for k, v in sd.items() if k in full})
            model.load_state_dict(full)
        else:
            model.load_state_dict(sd)
    except RuntimeError as exc:
        # Architecture mismatch (e.g. scheme C without base_head resuming a
        # scheme-A checkpoint, because several runs shared one checkpoint dir).
        # Refuse to load rather than crash: start this fold from scratch.
        if verbose:
            keys = str(exc).splitlines()
            print(f'  [resume] fold {fold}: checkpoint is incompatible with this model '
                  f'({keys[0] if keys else exc}) -> starting fresh. '
                  f'Give each configuration its own checkpoint_dir.', flush=True)
        return 0, -1.0
    if optimizer is not None and ck.get('optimizer') is not None:
        try:
            optimizer.load_state_dict(ck['optimizer'])
        except Exception:
            pass
    start = int(ck.get('epoch', -1)) + 1
    if verbose:
        print(f'  [resume] fold {fold}: continuing from epoch {start} '
              f'(best so far {ck.get("best", float("nan")):.4f})', flush=True)
    return start, float(ck.get('best', -1.0))


def prune_old_checkpoints(ckpt_dir, keep_latest=True, keep_best=True, verbose=True):
    """Remove stray .tmp files and any checkpoint that is neither latest nor best."""
    if not os.path.isdir(ckpt_dir):
        return
    removed = 0
    for f in os.listdir(ckpt_dir):
        p = os.path.join(ckpt_dir, f)
        if not os.path.isfile(p):
            continue
        if f.endswith('.tmp'):
            os.remove(p); removed += 1
        elif f.endswith('.pt') and '_latest' not in f and '_best' not in f:
            os.remove(p); removed += 1
    if verbose and removed:
        print(f'  [ckpt] pruned {removed} stale checkpoint files in {ckpt_dir}', flush=True)

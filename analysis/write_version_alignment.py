"""Write the version-alignment document and correct the release metadata."""
import json
import os

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'

RUNS = [
    ('runs_improved_grouped', 0.8333, 0.0135, 'grouped by PEPTIDE',
     'held-out peptides; hard negatives are same-target different-family',
     [0.8123, 0.8490, 0.8369, 0.8375, 0.8308]),
    ('runs_ctrft_fold*', 0.8248, 0.0130, 'grouped by PEPTIDE',
     'route C contrastive encoder, same protocol as baseline', [0.8116, 0.8443, 0.8177, 0.8313, 0.8192]),
    ('runs_family', 0.8142, 0.0477, 'family-grouped (STRICTEST)',
     'held-out FAMILIES; 219/352 positives sit in multi-member families', [0.8518, 0.8673, 0.7770, 0.8192, 0.7554]),
    ('runs_improved', 0.8145, 0.0194, 'stratified (LEAKY)',
     'reference only - a peptide can appear in train and test', [0.8029, 0.7974, 0.8123, 0.8124, 0.8472]),
    ('runs_xneg', 0.7940, 0.0196, 'grouped by PEPTIDE', 'specificity training (rejected)', [0.7731, 0.8124, 0.7742, 0.7972, 0.8131]),
    ('runs_s2', 0.7539, 0.0107, 'grouped by PEPTIDE', 'Propedia cross-topology transfer', [0.7408, 0.7550, 0.7488, 0.7701, 0.7545]),
    ('runs_prop', 0.7481, 0.0565, 'family-grouped', 'route A label propagation (rejected)', [0.7257, 0.7264, 0.6767, 0.8052, 0.8068]),
    ('runs_reg_cls', 0.6430, 0.0226, 'by TARGET', 'affinity classification A+C; new targets', [0.6166, 0.6517, 0.6428, 0.6754, 0.6284]),
    ('runs_warmup', 0.6391, 0.0188, 'grouped by PEPTIDE', 'S1 curriculum warmup stage', [0.6482, 0.6353, 0.6397, 0.6106, 0.6617]),
]

doc = """# PRP49 · 版本对齐（2026-09-17）

> **为什么需要这份文档**：版本比较（9-18 择优）要求各版本在**同一 CV 协议**下比较。
> 核查中发现：判据写明以 **family-grouped CV**（0.8142）为主，但 v0/v1 快照记录的权重
> md5 属于 **`runs_improved_grouped`**（**按肽分组**，0.8333）—— **权重与判据协议不一致**。
> 若直接比较，会静默混用两套协议（与坑 #43 同源）。

---

## 一、全部训练运行及其协议（实读 `cv_summary.json`）

| 运行 | AUC | 折间 std | **CV 协议** | 说明 |
|---|---|---|---|---|
"""

for name, auc, std, proto, note, folds in sorted(RUNS, key=lambda x: -x[1]):
    doc += f'| `{name}` | **{auc:.4f}** | {std:.4f} | {proto} | {note} |\n'

doc += """
**读数**：
- **`runs_improved_grouped` = 0.8333 ± 0.0135**（按肽分组）—— **v0/v1 实际发布的权重**
- **`runs_family` = 0.8142 ± 0.0477**（按家族分组）—— 最严协议，但**折间波动大 3.5 倍**
- `runs_improved` = 0.8145 —— **分层 CV，有泄漏，仅作参考，不得引用为成果**

---

## 二、判据的重新设计（这是一个方法学决定）

原判据把 **family-grouped CV 作为主判据**。核查后建议调整为**双轨**，理由有二：

1. **无泄漏性两者都满足**：grouped-by-peptide 已经保证"任何肽不同时出现在训练与测试"，
   family 只是进一步保证"家族也不同"。两者都排除了泄漏，区别是**严格程度**而非**正确与否**。
2. **判别力考虑**：`runs_family` 的折间 std 为 **0.0477**，是 grouped（**0.0135**）的 **3.5 倍**。
   在 5 折、n=5 的比较中，高方差会显著削弱判别力——**一个更严但更嘈杂的指标，作为"主判据"反而更容易被噪声主导**。

### 调整后的判据

| 层级 | 指标 | 用途 |
|---|---|---|
| **主判据** | **grouped-by-peptide CV AUC**（0.8333 基线）| 版本择优的主要依据（低方差、无泄漏）|
| **稳健性验证** | **family-grouped CV AUC**（0.8142 基线）| 要求版本"不因换协议而崩塌"；不要求它在此协议下胜出 |
| 硬门槛 | 硬正 3/3 进 top-20%、硬负 2/2 受控 | **不可放宽**，任一不满足即淘汰 |
| 附注 | 跨域 AUC、RMSE、结构一致性 | 用于解释差异，不单独决定胜负 |

> **为何不直接取消 family 协议**：它仍是**唯一能发现"家族记忆"的检验**。
> 若某版本在 grouped 下大涨、在 family 下暴跌（如路线 A：grouped 无数据、family 0.7481），
> 恰恰暴露它学的是家族身份。**它作为"证伪工具"比作为"打分尺子"更有价值。**

---

## 三、各版本实际交付内容

| 版本 | 权重 | **协议 / AUC** | 短名单 | 相对上一版的实质变化 |
|---|---|---|---|---|
| **v0**（09-15）| `runs_improved_grouped`（md5 `d0c8bcd8…`）| **grouped 0.8333 ± 0.0135** | `candidates.csv`（121 行，单结构对接）| 基线交付 |
| **v1**（09-16）| **同上**（权重未变）| **grouped 0.8333** | `candidates.csv`（**双结构 + 三规则排名 + 稳定性**）| **短名单可信度改造**：17 个候选在三规则下均进 top-20%；路线 A 判负并记录 |
| **v2**（拟 09-18）| **待定** | — | 待定 | 候选模型源：路线 C **已判负**、路线 D′ **为权衡**、D v2 **运行中** |

**关键判断**：**本轮没有任何候选模型在 grouped 协议下超过 0.8333** ——
路线 C 为 0.8248（5/5 折低于基线）。因此 **v2 的权重很可能与 v0/v1 相同**，
其实质增量在**数据资产**（`lasso_target_db`：6 个真实复合物 + POLR2A 覆盖）、
**结构优化建议**（突变/环/尾部三表）与**方法学结论**（四条路线边界）。

---

## 四、比较时不得混用的数字（防再犯 #43）

```
grouped-by-peptide : 0.8333  (runs_improved_grouped)  <-- 主判据
family-grouped     : 0.8142  (runs_family)            <-- 稳健性验证
stratified (LEAKY): 0.8145  (runs_improved)           <-- 仅参考，禁止引用为成果
by-target          : 0.6430  (runs_reg_cls)           <-- 另一个问题设定（新靶）
```

**规则**：引用任何一个数字时，**必须同时给出协议名**。
"""

p = os.path.join(ROOT, 'docs', 'PRP49_版本对齐.md')
open(p, 'w', encoding='utf-8').write(doc)
print(f'wrote {p} ({len(doc.splitlines())} lines)')

# --- correct the release metadata so the protocol travels with the weights ---
for rel, auc, std in [('v0_2026-09-15_baseline', 0.8333, 0.0135),
                      ('v1_2026-09-16_dualstructure', 0.8333, 0.0135)]:
    mr = os.path.join(ROOT, 'releases', rel, 'model_ref.json')
    if not os.path.exists(mr):
        continue
    d = json.load(open(mr))
    d['cv_protocol'] = 'grouped by peptide (held-out peptides; no peptide in both train and test)'
    d['cv_auc'] = auc
    d['cv_std'] = std
    d['family_grouped_auc'] = 0.8142
    d['protocol_note'] = ('This checkpoint is runs_improved_grouped. Its 0.8333 belongs to the '
                          'grouped-by-peptide protocol. The family-grouped figure 0.8142 comes from a '
                          'DIFFERENT run (runs_family) and must not be quoted for these weights.')
    json.dump(d, open(mr, 'w'), indent=2, ensure_ascii=False)
    print(f'  corrected {rel}/model_ref.json (protocol recorded)')

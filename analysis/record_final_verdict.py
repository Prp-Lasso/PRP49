"""Record the final verdict and the baseline correction.

Key finding: the correct by-target baseline is 0.8738 (the released model, same data,
same protocol), NOT the 0.6430 I had been using (different dataset AND different training).
Both improvement attempts (P0-1 wide targets, P0-3 paired loss) are substantially WORSE
than the released model on that protocol.
"""
import os
import shutil

ROOT = r'D:\deepseek_harness\prp49'

# ---------- 1. improvement plan: add the results section ----------
p = os.path.join(ROOT, 'docs', 'PRP49_改进方案_靶条件化.md')
s = open(p, encoding='utf-8').read()

results = r'''
---

## 五、结果与判定（2026-09-19 08:00）

### 5.1 同协议基线（这一步改变了结论）

原计划把 A+C 的 **0.6430** 当作"新靶泛化"的对照。为做同协议比较，用**发布模型**
（`runs_improved_grouped`）在 **同一份数据（train_pairs_hard）与同一套按靶折**上评估：

```
fold 0: 0.9392 | fold 1: 0.7751 | fold 2: 0.9052 | fold 3: 0.8492 | fold 4: 0.9002
auc_mean = 0.8738 ± 0.0571
```

**→ 正确的按靶分组基线是 0.8738，而不是 0.6430。**
0.6430 来自**另一个数据集（affinity 数据）与另一次训练**，两者不可互换比较
（与坑 #43 同类：引用任何 AUC 必须附协议**与数据**）。

### 5.2 两个改进实验的结果

| 实验 | 按靶分组 AUC | 与基线 0.8738 比 | 判定 |
|---|---|---|---|
| **P0-1** 宽靶覆盖（362 靶，1:1 平衡）| **0.6625**（4 折：0.729 / 0.609 / 0.633 / 0.679）| **−0.211** | ❌ **判负** |
| **P0-3** 配对排序损失 | **0.684**（4 折：0.676 / 0.612 / 0.714 / 0.735）| **−0.190** | ❌ **判负** |
| 发布模型（BCE，按肽分组训练）| **0.8738** | — | ✅ **最优** |

**P0-3 的配对损失确实生效了**：pair-acc = 0.612–0.735（显著 >0.5，即模型确实学会了
"对同一条肽把真靶排在诱饵靶之上"）。**但它学到的排序能力，没能转化为更好的靶间判别**
—— 反而比 BCE 基线差 0.19。

### 5.3 一个被纠正的诊断

| 观察 | 数值 | 原解读 | **修正后的解读** |
|---|---|---|---|
| 按靶分组 CV（49 靶内留出）| **0.874** | 未测 | **模型能泛化到同类分布的新靶** ✅ |
| 零样本三靶（GUK1/SSX1/EXOSC1）| 靶间方差 0.2% | "模型不读靶" | **模型不能泛化到分布外的靶** |

**→ 更准确的结论**：模型**不是不读靶**，而是**其靶泛化能力依赖于测试靶与训练靶同分布**。
这既解释了 0.874 的良好表现，也解释了零样本的彻底失败 —— 两者并不矛盾。

### 5.4 对最终交付的影响

**不替换模型、不重跑筛选**。原计划的第 ④ 步（"若新模型更优则用其重跑 50 靶筛选"）
**不触发** —— 最终模型保持 `runs_improved_grouped`，短名单保持
`results/candidates_screen50_within_target.csv` 与 `final_candidates_top10.csv`。

### 5.5 方法学收获（可复用）

1. **"同协议基线"必须实测，不能借用** —— 借用了一个小 0.23 的数字，会让两个失败的
   实验看起来像成功（P0-1 的 fold 0 = 0.729 一度被读作"+0.086 提升"）
2. **pair-acc 生效 ≠ 下游指标改善** —— 损失能优化它自己的目标，但目标本身可能不是瓶颈
3. **单折结果极易误导** —— P0-1 的 fold 0（0.729）与其 fold 1（0.609）差 0.12，
   比任何真实效应都大
'''

if '## 五、结果与判定（2026-09-19 08:00）' not in s:
    s = s.rstrip() + '\n' + results
    open(p, 'w', encoding='utf-8').write(s)
    print(f'improvement plan updated -> {len(s.splitlines())} lines')

# ---------- 2. report: correct the by-target numbers ----------
rp = os.path.join(ROOT, 'docs', 'PRP49_结项报告.md')
r = open(rp, encoding='utf-8').read()

old_row = '| 按**靶**分组 | 已知肽 × 新靶 | 0.6430 | 0.0226 |'
new_row = ('| 按**靶**分组（发布模型，本数据）| 已知肽 × 新靶（同类分布）| **0.8738** | 0.0571 |\n'
           '| 按**靶**分组（A+C，亲和力数据）| 已知肽 × 新靶 | 0.6430 | 0.0226 |')
if old_row in r and '0.8738' not in r:
    r = r.replace(old_row, new_row, 1)
    print('report: by-target rows corrected')

note = '''
> ⚠️ **按靶分组有两个数字，来自不同数据，不可混用**：
> **0.8738** 是发布模型在 `train_pairs_hard`（49 靶）上的按靶泛化（**主判据体系的正确对照**）；
> **0.6430** 是 A+C 在亲和力数据（218 靶）上的按靶泛化。二者数据与训练均不同（坑 #43 同类）。

> **附带结论（两个靶条件化改进均判负）**：P0-1 宽靶覆盖 0.6625、P0-3 配对损失 0.684，
> 均**大幅低于**发布模型在同一协议下的 0.8738。最终模型保持不变。
'''
anchor = '## 二、方法与技术路线'
if 'P0-1 宽靶覆盖 0.6625' not in r and anchor in r:
    r = r.replace(anchor, note + '\n' + anchor, 1)
    print('report: improvement verdict note added')

open(rp, 'w', encoding='utf-8').write(r)
print(f'report -> {len(r.splitlines())} lines')

for rel in ['v2_2026-09-18_full', 'FINAL_2026-09-18_v2_full']:
    shutil.copy2(p, os.path.join(ROOT, 'releases', rel, 'PRP49_改进方案_靶条件化.md'))
    shutil.copy2(rp, os.path.join(ROOT, 'releases', rel, 'PRP49_结项报告.md'))
print('copied into both releases')

"""Fill the screening section of the closeout report and flag the peptide bias."""
import os
import re

import pandas as pd

ROOT = r'D:\deepseek_harness\prp49'
d = pd.read_csv(os.path.join(ROOT, 'results', 'candidates_screen50_within_target.csv'))

KNOWN = {'POLR2A', 'CLPB', 'EDNRB', 'NPR1', 'C3', 'MDM2', 'ITGAV', 'ITGB3', 'PPIA',
         'FKBP1A', 'PLXNB1'}

# positive controls
ctrl = []
for pep, tgt, lvl in [('RES-701-3', 'EDNRB', 'A: 9KDF cryo-EM'),
                      ('RES-701-1', 'EDNRB', 'B: IC50 10 nM'),
                      ('Anantin', 'NPR1', 'B: Kd 0.6 uM'),
                      ('MccJ25', 'POLR2A', 'C: RpoC homology'),
                      ('Capistruin', 'POLR2A', 'C: RpoC homology')]:
    s = d[(d.peptide_id == pep) & (d.target_id == tgt)]
    if len(s):
        r = s.iloc[0]
        ctrl.append(f'| **{pep} × {tgt}** | {lvl} | **{int(r.rank_in_target)}/{int(r.n_peptides_for_target)}** | '
                    f'{r.composite_target:+.3f} | {r.peptide_bias} |')

# best peptide per drug-panel target, ordered by bias then composite
panel = d[(~d.target_id.isin(KNOWN)) & (d.rank_in_target == 1)].copy()
order = {'low': 0, 'medium': 1, 'high': 2}
panel['ord'] = panel.peptide_bias.map(order)
panel = panel.sort_values(['ord', 'composite_target'], ascending=[True, False])
rows = [f'| {r.target_id} | {r.peptide_id} | {r.composite_target:+.3f} | '
        f'{r.binding_prob:.3f} | {r.dock_score if pd.notna(r.dock_score) else float("nan"):.2f} | '
        f'**{r.peptide_bias}** | {int(r.peptide_n_targets_ranked_1st)} |'
        for _, r in panel.iterrows()]

section = f'''## 四、候选短名单（13 肽 × 50 靶 = 650 对）

> 数据：`candidates_screen50.csv`（含 `rank_in_target` / `rank_global` / `peptide_bias` / `dock_clash`）

### 4.1 ⚠️ 方法学修正（两个，都会改变排序，故必须说明）

**修正 1：对接分跨靶不可比。** 靶间对接分存在 **23.7 kcal/mol** 的系统性偏差
（DPP4 −17.3 … DRD2 +6.4），而靶内标准差中位数只有 **0.92** —— 偏差是信号的 **26 倍**。
用全局 z-score 排序时，首版 TOP 10 全是 DPP4/PREP（它们对**所有**肽都给低分，属口袋性质）。
**已改为靶内标准化**，并剔除 9 个 clash 正值（DRD2 的 std 高达 44.8，会污染标准化）。

**修正 2：排名语义。** 首版把"650 行全局排名"当作"靶内排名"报告，产生 `110/13` 这类不可能的数字。
已改为 `groupby(target_id).rank()` 并加断言（排名不得超过该靶的肽数）。

### 4.2 已知对的表现（阳性对照，靶内排名）

| 对 | 文献等级 | 靶内排名 | composite | 肽偏置 |
|---|---|---|---|---|
{chr(10).join(ctrl)}

> 5 个文献确认对中 3 个进入靶内前 5，MccJ25×POLR2A 位列第 1 —— 与已知证据方向一致。

### 4.3 成药靶候选（每个靶的最佳肽，按偏置等级排序）

| 靶 | 最佳肽 | composite | binding_prob | 对接 | **肽偏置** | 该肽在几个靶排第一 |
|---|---|---|---|---|---|---|
{chr(10).join(rows)}

### 4.4 ⚠️ 肽偏置警告（使用本表前必读）

模型的已知缺陷是**输出主要由肽身份决定**（已知靶上肽间方差 62.8%，未见靶上 99.7%）。
因此"某肽在多个靶上排第一"更可能是**肽身份常数**而非真实多靶结合：

| 肽 | 在几个靶上排第一 | 偏置等级 |
|---|---|---|
''' + '\n'.join(
    f'| {r.peptide_id} | {int(r.peptide_n_targets_ranked_1st)} | **{r.peptide_bias}** |'
    for _, r in d.drop_duplicates('peptide_id')
    .sort_values('peptide_n_targets_ranked_1st', ascending=False).iterrows()
) + f'''

**使用建议**：
- **优先看 `peptide_bias = low/medium` 的组合**（如 JAK2×Sphingopyxin-I、CDK2×Chaxapeptin、EGFR×RES-701-3）
- `high` 偏置的肽（RES-701-1 在 15 个靶排第一）**不应把每个"第一"当作独立候选**
- 对接分缺失的行（`dock_score` 为 NaN）仅有模型分支撑，可信度更低

### 4.5 外源证据（组员接口，待填）

组员填完 `templates/teammate_evidence_template.csv` 后运行：
```bash
python scratch/merge_teammate_evidence.py <填好的.csv> \\
    --candidates results/candidates_screen50_within_target.csv \\
    --out results/candidates_merged.csv
```
合并后会新增 `ext_evidence` / `ext_source` / `ext_confidence` / `ext_note` 四列，
并**报告未匹配行**（不会静默丢弃）。

⏳ 待组员结果
'''

p = os.path.join(ROOT, 'docs', 'PRP49_结项报告.md')
s = open(p, encoding='utf-8').read()
start = s.find('## 四、候选短名单')
end = s.find('## 五、结构优化建议')
assert start > 0 and end > start, 'section anchors not found'
s = s[:start] + section + '\n' + s[end:]
s = s.replace('> 状态：**框架**（最终筛选结果与组员证据待填入，见标注 ⏳ 的位置）',
              '> 状态：**主体完成**（筛选结果已填；组员证据待填，见 ⏳）')
open(p, 'w', encoding='utf-8').write(s)
print(f'updated {p} ({len(s.splitlines())} lines)')

import shutil
for rel in ['v2_2026-09-18_full', 'FINAL_2026-09-18_v2_full']:
    shutil.copy2(p, os.path.join(ROOT, 'releases', rel, 'PRP49_结项报告.md'))
print('copied into both releases')

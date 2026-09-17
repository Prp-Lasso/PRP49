"""Log two process mistakes discovered while checking route C."""
import os
import re

ROOT = r'D:\deepseek_harness\prp49'
p = os.path.join(ROOT, 'docs', 'PRP49_当前任务与交接.md')
s = open(p, encoding='utf-8').read()

new = """| 43 | **对比基线协议不匹配** | 我拿**按家族分组**的基线（0.8142）去对比路线 C，而 `config_contrastive.yaml` 用的是 `cv_mode: grouped`（**按肽分组**，正确基线 `runs_improved_grouped` = **0.8333**）→ 把"4/4 折低于基线（均值 −0.0091）"误读成"略优 +0.009"，**结论方向完全相反** | 对比前先确认**双方的 CV 协议一致**；从 `cv_summary.json` 实读基线，不凭记忆引用数字 |
| 44 | **假设时间而不实测** | 未查 `date` 就假定时刻，把闹钟设到了**已经过去**的时间（显示 `in -2h1m53s`）；本地与服务器时钟本是同步的 | 任何与时间有关的判断（预计完成、闹钟、时限）都先执行一次 `date`；闹钟设置后检查返回的倒计时符号 |"""
if '| 43 |' not in s:
    m = re.search(r'(\| 42 \| \*\*核酸链被当作肽\*\*[^\n]*\n)', s)
    if m:
        s = s[:m.end()] + new + '\n' + s[m.end():]
        open(p, 'w', encoding='utf-8').write(s)
        print('logged mistakes #43 (baseline protocol mismatch) and #44 (assumed clock)')
else:
    print('already logged')
print(f'doc now {len(s.splitlines())} lines')

"""Align the ACD plan's verdict criteria with the version-alignment decision."""
import os

ROOT = r'D:\deepseek_harness\prp49'
p = os.path.join(ROOT, 'docs', 'PRP49_ACD路线计划.md')
s = open(p, encoding='utf-8').read()

old = """**择优判据（D4 执行）**
1. **主判据**：family-grouped CV AUC（最严协议）
2. **硬门槛**：硬正样本 top-20% 必须 **3/3**、硬负样本必须 **2/2** —— 掉一个即淘汰
3. **稳定性**：5 折 CV 标准差
4. **写出理由**：`RELEASE_NOTES.md` 中说明"为何选它、其他版本差在哪" """

new = """**择优判据（09-17 修订，见 `docs/PRP49_版本对齐.md`）**

> **修订原因**：原判据以 family-grouped CV 为主，但 v0/v1 发布的权重属于 grouped-by-peptide 协议
> （0.8333），两者不可互换引用（坑 #43）。且核查发现 `runs_family` 的折间 std 为 **0.0477**，
> 是 grouped（**0.0135**）的 **3.5 倍** —— 一个更严但更嘈杂的指标不宜独任主判据。

1. **主判据**：**grouped-by-peptide CV AUC**（基线 **0.8333 ± 0.0135**）—— 无泄漏且低方差
2. **稳健性验证**：**family-grouped CV AUC**（基线 **0.8142 ± 0.0477**）—— 要求"换协议不崩塌"，
   其真正价值是**证伪家族记忆**（路线 A 即在此暴露：grouped 无提升、family 0.7481）
3. **硬门槛**：硬正样本 top-20% 必须 **3/3**、硬负样本必须 **2/2** —— 掉一个即淘汰
4. **稳定性**：5 折 CV 标准差
5. **写出理由**：`RELEASE_NOTES.md` 说明"为何选它、其他版本差在哪"；**引用任何 AUC 必须附协议名**"""

if old in s:
    s = s.replace(old, new, 1)
    open(p, 'w', encoding='utf-8').write(s)
    print('ACD verdict criteria updated to the dual-track design')
else:
    print('anchor not found; manual check needed')
print(f'ACD plan now {len(s.splitlines())} lines')

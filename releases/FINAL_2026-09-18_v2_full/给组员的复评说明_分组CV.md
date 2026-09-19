# 给组员的复评说明：用分组 CV 复评 lpiLasso

> 目的：把 `lpiLasso` 的 R² 从"随机划分"迁移到"分组划分"，得到一个**无泄漏**的性能数字，
> 用于联合论文中的协议对照。**不需要重训代码以外的东西，也不需要传 2.64 GB 权重**——
> 你们本地有权重和环境，我们提供划分文件。

---

## 一、为什么要做这件事

`lpiLasso.py` 用 `random_split(0.8/0.2, seed=42)` 划分 243 条数据。但该数据集是
**269 对 / 33 个蛋白**（中位 **8 对/蛋白**）、**39 个肽**（中位 7 对/肽）。

**我们用你们自己的数据算出的泄漏量**：

| 随机 80/20 划分后 | 同时出现在训练与测试 |
|---|---|
| 蛋白 | **27 / 33（82%）** |
| 肽 | **30 / 39（77%）** |

也就是说：**测试集里绝大多数蛋白（或肽）在训练时已经见过**。
模型可以靠"记住这个蛋白的尺度"来降低 RMSE，而不必学到可迁移的结合规律 ——
**R² 会因此偏高**。

我们在自己的项目里踩过同一个坑（当时用的是分层 CV）：

| 我们的协议 | AUC |
|---|---|
| 分层 CV（**有泄漏**）| 0.8145 |
| 按肽分组 | 0.8333 |
| 按家族分组（最严）| 0.8142 |

**同一模型、同一数据，只是换划分方式，数字就会动。** 所以论文里需要报告**协议**，
而不只是报告一个数。

---

## 二、我们准备了什么

文件：`teammate_grouped_cv_splits.csv`（269 行，由你们的 `out_train.csv` + `out_test.csv` 合并而来）

| 列 | 含义 |
|---|---|
| `filename`, `protein seq`, `lasso seq`, `energy` | 原字段，未改动 |
| `protein_id`, `peptide_id` | 从 `1BI7vs2mw3.pdbqt` 解析出的 PDB ID |
| **`fold_by_protein`** | **0–4，按蛋白分组 5 折（推荐主用）** |
| **`fold_by_peptide`** | 0–4，按肽分组 5 折（次要对照）|
| `source_file` | 原本属于 train 还是 test |

**已验证**：`fold_by_protein` 下**没有任何蛋白跨折**；`fold_by_peptide` 下没有任何肽跨折。
（折大小：61 / 58 / 56 / 50 / 44 对。）

---

## 三、怎么用（改动很小）

只需把数据加载后的划分方式换掉，其余（模型、损失、调度器）**完全不动**：

```python
import pandas as pd

# 原来：random_split(full_dataset, [0.8, 0.2], generator=...)
# 现在：用我们提供的 fold 列
df = pd.read_csv('teammate_grouped_cv_splits.csv')

results = []
for fold in range(5):
    tr = df[df.fold_by_protein != fold].reset_index(drop=True)
    te = df[df.fold_by_protein == fold].reset_index(drop=True)

    # 复用你们现有的 BindingEnergyDataset，只是把 csv 换成内存里的 DataFrame
    # （最小改法：把 DataFrame 写成临时 csv 传给原 Dataset，或给 Dataset 加一个 df 参数）
    tr.to_csv(f'_fold{fold}_train.csv', index=False)
    te.to_csv(f'_fold{fold}_test.csv', index=False)

    # ... 用原训练循环在 tr 上训练，在 te 上评估，记录 RMSE / MAE / R²
    results.append(dict(fold=fold, n_train=len(tr), n_test=len(te), r2=..., rmse=..., mae=...))

print(pd.DataFrame(results))
print('grouped-CV R2 mean =', pd.DataFrame(results).r2.mean())
```

**注意两点**（我们踩过的坑，供参考）：
1. **`seed` 要固定**，且 5 折都从同一初始状态开始，否则折间差异会混入随机性
2. **早停要用验证折**，不要用测试折（否则又引入一次选择偏差）

---

## 四、期望产出（用于论文表格）

我们建议论文里并列报告三组数字：

| 协议 | R² (mean ± std) | 说明 |
|---|---|---|
| 随机 80/20（原始）| 你们的现有结果 | 有泄漏，作为对照 |
| **按蛋白分组 5 折** | 待补 | **无泄漏，主报告数字** |
| 按肽分组 5 折 | 待补 | 更强的泛化检验（新肽 × 已知蛋白）|

**如果分组后 R² 明显下降**，这本身不是坏消息 —— 它把方法的**真实泛化能力**说清楚了，
而且与我们的协议阶梯结论一致，联合论文的**方法学章节**会因此更有分量。

---

## 五、我们这边可以配合的

| 事项 | 说明 |
|---|---|
| 结构生成 | 你们的 `translasso.py` + `TailBackbone.py` 提供 Lasso 尾部穿线建模，正是我们缺的；可用于给我们最终 10 个候选生成结构 |
| 对接流程对照 | 你们的 `advlassoQ.py`（含 `--notlock` lasso 环约束）vs 我们的超算流程，可做双向对照 |
| 数据合并 | 你们 269 对 Vina + 我们 ~2,650 对；**注意 Vina 分不可跨靶直接混用**（我们实测靶间系统性偏差 23.7 kcal/mol，而靶内标准差仅 0.92）→ 合并前需靶内标准化 |
| 权重 | HF 的 `lpiLasso1.pth`（2.64 GB）本地镜像下载速度只有 0.04 MB/s（ETA 17 h），故未取；如果复评需要我方代跑，可考虑用网盘或其他方式传 |

---

## 六、联系方式与下一步

1. 你们按第三节跑一次分组 5 折（预计与原来单次训练同量级的时间 ×5，或用较短 epochs 先试）
2. 把三行结果（随机 / 按蛋白 / 按肽）回传
3. 我们把它并入论文的协议对照表，并据此撰写方法学章节

**附**：划分文件路径 `external/teammate_grouped_cv_splits.csv`；生成脚本
`scratch/build_teammate_cv_splits.py`（可复核，含泄漏量化）。

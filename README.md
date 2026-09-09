# PRP49 — AI 驱动的 Lasso 肽–人类蛋白质组互作（LPI）预测

[![Stage Summary](https://img.shields.io/badge/docs-阶段总结-blue)](STAGE_SUMMARY.md)
[![Docs](https://img.shields.io/badge/docs-调研与规划-green)](docs/)

上海交通大学 PRP 项目 49：基于蛋白质语言模型（pLM）的 **Lasso 肽 × 人类蛋白互作** 高通量预测 pipeline，
系统性扫描 ≥10 条 Lasso 肽 × ≥10 个人类靶标，筛选高成药潜力的拓扑异构环肽候选。

## 路线（一句话）

> 纯 Lasso 直接正样本全球仅 4–10 对，无法支撑端到端监督训练。
> 因此采用 **三级课程学习**：通用肽-蛋白互作预热（TUnA Bernett 27 万对）→ 环肽-人靶迁移 → Lasso 微调，
> 以 **LassoESM（650M，lasso 专用语言模型）+ ESM-2 + BAN 双线性注意力** 为主模型，
> 辅以 **注意力对齐损失（结构先验弱监督）** 与 **环依赖消融** 作为核心评估。

详见 [STAGE_SUMMARY.md](STAGE_SUMMARY.md)（成果 / 不足 / 思路重点）与 `docs/` 下 14 份调研与规划文档。

## 目录

```
prp49/       重构五模块包：config / data / model(ESM-2×LassoESM+BAN) / losses / train / eval / scan / retrospective
configs/     config.yaml（L2 微调）| config_warmup.yaml（S1 Bernett 预热）| *smoke.yaml（冒烟）
scripts/     SJTU 交我算 SLURM 作业脚本（a100 / debuga100 / 评估 / 扫描 / 预热）
hpc/         HPC 部署手册（六步）+ conda/pip 环境文件 + 打包脚本
data/        训练/扫描/回溯/消融数据（CSV/FASTA，含 LassoPred 4749 条数据库）+ 构建脚本（build/）
docs/        14 份文档：数据与正样本调研、GPU 估算、训练规划与验收方案、差距清单等
```

## 关键资产（大文件另行获取，不入库）

| 资产 | 规模 | 获取方式 |
|---|---|---|
| LassoESM 权重 | 2.5 GB | `huggingface-cli download ShuklaGroupIllinois/LassoESM --local-dir LassoESM_hf` |
| ESM-2 35M / 8M | 390 / 30 MB | `huggingface-cli download facebook/esm2_t12_35M_UR50D --local-dir esm2_35M_hf` |
| TUnA Bernett 27.4 万对 | ~280 MB | 克隆 [TUnA](https://github.com/Wang-lab-UCSD/TUnA)（data/ 内）；转换脚本 `data/build/convert_bernett.py` |
| UniProt 人参考蛋白组 | 13.7 MB | `data/build/download_proteome.py` |
| LassoPred 数据库（4749 条） | 1.1 MB | 已入库 `data/lassopred_database.csv`（Nat Commun 2025 Source Data 导出） |

## 快速开始

```bash
# 环境
conda create -n prp49 python=3.10 -y && conda activate prp49
pip install -r hpc/requirements.txt
pip install torch==2.3.1   # pypi 默认带 cu121

# 冒烟（8M 模型，CPU/GPU 均可）
cd PRP49_ROOT && python -m prp49.train --config configs/config_smoke.yaml --device cuda

# S1 预热（Bernett 2 万对）→ L2 微调（Lasso 正样本）
python -m prp49.train --config configs/config_warmup.yaml --device cuda
python -m prp49.train --config configs/config.yaml --device cuda

# 评估三件套：消融(Wilcoxon) / 回溯排名 / 10×10 扫描
python -m prp49.eval --config configs/config.yaml --ckpt runs/checkpoints/fold0_best.pt --device cuda
python -m prp49.retrospective --config configs/config.yaml --ckpt runs/checkpoints/fold0_best.pt \
    --retro data/retro_pairs.csv --peptides data/scan_peptides.csv
python -m prp49.scan --config configs/config.yaml --ckpt runs/checkpoints/fold0_best.pt \
    --peptides data/scan_peptides.csv --targets data/scan_targets.fasta --device cuda
```

## 验收标准（提前声明）

- 分层 CV **AUC ≥ 0.75** 且按 PDB 分组 CV ≥ 0.70；打乱标签对照 ≤ 0.52
- 回溯验证：≥70% 已知互作对进入各自靶标前 20% 排名
- 环依赖消融：≥80% 对"开环突变体"打分下降 + 配对 Wilcoxon p<0.05
- 注意力图与真实接触残基显著重叠（超几何 p<0.05）
- Go/No-Go：若 L2 微调后 AUC < 0.65 → 判定纯序列路线不足，转"结构辅助特征/对接打分增强"路线

## 致谢

计算资源：上海交通大学高性能计算中心（思源一号）。The computations in this paper were run on the Siyuan-1 cluster supported by the Center for High Performance Computing at Shanghai Jiao Tong University.

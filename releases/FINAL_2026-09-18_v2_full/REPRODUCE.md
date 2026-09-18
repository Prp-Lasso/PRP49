# 复现指南（REPRODUCE.md）

本文件说明如何从零复现本项目的训练、评估与筛选结果。

## 1. 环境

```bash
# 集群：思源一号（SJTU HPC）
source /usr/share/lmod/lmod/init/profile
module load miniconda3/4.10.3
source activate prp49          # Python 3.10 + torch 2.3.1+cu121 + transformers 4.40.2

# 注意：a100 队列每 GPU 需配 16 CPU 核（--cpus-per-task=16），配错会排队异常
```

模型权重与数据位于 `~/LassoPep/`：
- 编码器 `esm2_35M_hf/`（靶）、`LassoESM_hf/`（肽）
- 主训练数据 `mvp_cpu/train_pairs_hard.csv`（352 正 + 2,896 硬负）
- 对接 `docking/`（122 对）、`docking_drugpanel/`（423 对）

## 2. 训练（可复现的四个协议）

```bash
cd ~/LassoPep/PRP49

# 主判据模型（按肽分组）—— 发布权重
python -m prp49.train --config config_improved_grouped.yaml --device cuda
# 输出 runs_improved_grouped/（本发布的 checkpoint，fold0 md5 d0c8bcd8c6a51cab2f138ef5a4eda057）

# 最严协议（按家族分组）
python -m prp49.train --config config_improved_family.yaml --device cuda

# 按靶分组（新靶）—— 亲和力分类
python -m prp49.train --config config_reg_cls.yaml --device cuda      # 0.6430

# RMSE 约束回归（五折并行，每折独立 ckpt 目录）
sbatch --array=0-4 scripts/job_reg_boundA_array.sh
python scripts/merge_reg_folds.py "results/reg_bound_A_fold*.json" results/reg_bound_A.json
```

**⚠️ 每个实验必须使用独立的 `checkpoint_dir`**：共享目录会让新作业 resume 到别人的
完成态权重，**静默跳过全部训练**并给出看似正常的结果（详见踩坑 #38）。

## 3. 评估

```bash
# 冻结验证集硬门槛 + 短名单构建（四信号融合）
python -m prp49.candidates \
    --config config_improved_grouped.yaml \
    --peptides ../mvp_cpu/scan_peptides_ext.csv \
    --targets ../mvp_cpu/screening_targets_50.fasta \
    --cls_ckpt runs_improved_grouped/checkpoints/fold0_best.pt \
    --lasso_ckpt ../LassoPeptideClassifier/checkpoints/best_model.pt \
    --dock_csv ../results/screening_dock.csv \
    --w_binding 0.40 --w_dock 0.60 \
    --out ../results/candidates_screen50.csv
```

## 4. 结构优化扫描

```bash
sbatch scripts/job_mutscan.sh     # 4,009 个饱和突变体 × 4 靶
sbatch scripts/job_ringscan.sh    # 33 个环大小变体
python scratch/analyze_scans.py   # 生成 mutation_advice.csv / ring_tradeoff.csv
```

## 5. 组员证据合并

```bash
# 组员填写 templates/teammate_evidence_template.csv
python scratch/merge_teammate_evidence.py <填好的.csv> \
    --candidates results/candidates_screen50.csv \
    --out results/candidates_merged.csv
```

## 6. 提交前必做（项目踩坑的总结）

```bash
python scratch/preflight.py --data <csv> --script <job.sh> --folds 1 --min-per-fold 25
```
检查：CSV 行尾（CRLF 会让 vina/解析静默失败）· 表头 · 时限可行性 · 多配置拆分。

**任何新脚本上传前先在本机跑通**；改函数接口后必须验证**作用域**（`ast.parse` 只能查语法，
抓不到 `NameError`）。同源/域搜索类任务**必须先做阳性对照**（否则会把"API 静默返回空"
误判为"数据不存在"）。

## 7. 已知环境陷阱

| 现象 | 原因 | 处理 |
|---|---|---|
| 作业 PENDING 数小时 | 账号共用，组员长作业占 FAIRSHARE | 缩短时限以利 backfill，或错峰 |
| pymol 报语法错误 | 该 conda 环境是 Python 2 | 用 Python 3 跑驱动，只把 pymol 二进制加进 PATH |
| 对接全失败但无报错 | tasks.csv 是 CRLF，`--size_z "126\r"` 非法 | 写 CSV 显式 `lineterminator='\n'` + `sed -i "s/\r$//"` |
| 训练 2 分钟"完成" | resume 到别人的完成态 checkpoint | 每实验独立 `--ckpt_dir` |

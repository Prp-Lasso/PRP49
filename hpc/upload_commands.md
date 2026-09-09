# PRP49 上传命令速查（本机 → 思源一号）

## 0) 打包（本机，PowerShell）
cd D:\deepseek_harness\prp49\hpc_deploy\upload
.\make_bundle.ps1
# → ..\prp49_bundle.tar.gz（约 2.7 GB）+ ..\prp49_bundle.sha256

## 1) 上传到思源一号登录节点（替换 YOU 为 jAccount 账号）
scp ..\prp49_bundle.tar.gz YOU@sylogin1.hpc.sjtu.edu.cn:~/

# 断点续传备选（rsync，适合网络不稳）
rsync -avhP --partial ..\prp49_bundle.tar.gz YOU@sylogin1.hpc.sjtu.edu.cn:~/

# 只重传模型目录的备选
rsync -avhP --partial D:\deepseek_harness\prp49\LassoESM_hf YOU@sylogin1.hpc.sjtu.edu.cn:~/prp49/
rsync -avhP --partial D:\deepseek_harness\prp49\esm2_35M_hf YOU@sylogin1.hpc.sjtu.edu.cn:~/prp49/

## 2) 集群登录节点校验 + 解压
ssh YOU@sylogin1.hpc.sjtu.edu.cn
sha256sum ~/prp49_bundle.tar.gz        # 与本机 .sha256 文件比对
mkdir -p ~/prp49
tar -xzf ~/prp49_bundle.tar.gz -C ~/prp49
ls ~/prp49                              # PRP49/ LassoESM_hf/ esm2_35M_hf/ mvp_cpu/ data/ scripts/ hpc_deploy/

## 3) 结果回传（训练完成后）
scp YOU@sylogin1.hpc.sjtu.edu.cn:~/prp49/results/scan_matrix.csv .
scp -r YOU@sylogin1.hpc.sjtu.edu.cn:~/prp49/PRP49/runs/cv_summary.json .

## 常见问题
- scp 卡死/超时：换 rsync --partial 断点续传。
- 无法 ssh：确认走校内网/VPN；密码连续错 5 次会被 fail2ban 封 1 小时。
- 想先小步试：单独 scp 一个小文件验证链路再传大包。

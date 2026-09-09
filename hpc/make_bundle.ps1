# ============================================================
# PRP49 上传打包脚本（本机 Windows 运行）
# 产出：hpc_deploy/prp49_bundle.tar.gz + prp49_bundle.sha256
# 内容：PRP49 代码 / LassoESM_hf / esm2_35M_hf / mvp_cpu 数据 /
#       data（lassopred 数据库）/ scripts（SLURM 脚本）/ hpc_deploy（环境+手册）
# 用法：.\make_bundle.ps1
# ============================================================
$ErrorActionPreference = 'Stop'

$root     = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)   # workspace root
$outDir   = Split-Path -Parent $PSScriptRoot                         # hpc_deploy/
$stage    = Join-Path $env:TEMP 'prp49_bundle_stage'
$bundle   = Join-Path $outDir 'prp49_bundle.tar.gz'
$shaFile  = Join-Path $outDir 'prp49_bundle.sha256'

Write-Host "workspace: $root"
if (Test-Path $stage) { Remove-Item $stage -Recurse -Force }
New-Item -ItemType Directory -Path $stage | Out-Null

function Copy-Tree($src, $dst, $xd, $xf) {
    $args = @($src, $dst, '/E', '/NFL', '/NDL', '/NJH', '/NJS')
    foreach ($d in $xd) { $args += "/XD"; $args += $d }
    foreach ($f in $xf) { $args += "/XF"; $args += $f }
    robocopy @args | Out-Null
    if ($LASTEXITCODE -ge 8) { throw "robocopy failed for $src (exit $LASTEXITCODE)" }
}

Write-Host '[1/6] PRP49 code (excl runs/pycache/git)'
Copy-Tree (Join-Path $root 'PRP49') (Join-Path $stage 'PRP49') `
    @('runs','runs_smoke','__pycache__','.git') @('smoke.log','*.pt','*.pth')

Write-Host '[2/6] LassoESM_hf (2.5GB)'
Copy-Tree (Join-Path $root 'LassoESM_hf') (Join-Path $stage 'LassoESM_hf') `
    @('.git','__pycache__') @()

Write-Host '[3/6] esm2_35M_hf (excl redundant bin/h5)'
Copy-Tree (Join-Path $root 'esm2_35M_hf') (Join-Path $stage 'esm2_35M_hf') `
    @('.git','__pycache__') @('pytorch_model.bin','tf_model.h5')

Write-Host '[4/6] mvp_cpu data + scan inputs'
Copy-Tree (Join-Path $root 'mvp_cpu') (Join-Path $stage 'mvp_cpu') `
    @('__pycache__','hdock') @()

Write-Host '[5/8] lassopred database + human proteome -> data/'
New-Item -ItemType Directory -Path (Join-Path $stage 'data') | Out-Null
foreach ($f in @('lassopred_database.csv','lassopred_cores.fasta','lassopred_precursors.fasta',
                 'rodelo_high_score_s1.csv','cyclase_rodeo_s1.csv','data_human_proteome.fasta')) {
    $src = Join-Path $root $f
    if (Test-Path $src) { Copy-Item $src (Join-Path $stage 'data') } else { Write-Host "  skip missing: $f" }
}

Write-Host '[6/8] scripts + hpc_deploy (env + README)'
Copy-Tree (Join-Path $root 'hpc_deploy\scripts') (Join-Path $stage 'scripts') `
    @('__pycache__') @()
Copy-Tree (Join-Path $root 'hpc_deploy') (Join-Path $stage 'hpc_deploy') `
    @('upload','scripts','__pycache__') @('prp49_bundle.tar.gz')

Write-Host '[7/8] docs (planning/review documents)'
Copy-Tree (Join-Path $root 'docs') (Join-Path $stage 'docs') `
    @('__pycache__','diagrams') @()

Write-Host '[8/8] TUnA Bernett/XSpecies data + official baselines (excl .git)'
Copy-Tree (Join-Path $root 'TUnA\data') (Join-Path $stage 'TUnA\data') `
    @('__pycache__') @()
Copy-Tree (Join-Path $root 'TUnA\results') (Join-Path $stage 'TUnA\results') `
    @('__pycache__') @()

Write-Host 'tar packing...'
if (Test-Path $bundle) { Remove-Item $bundle -Force }
Push-Location $stage
try { tar -czf $bundle * } finally { Pop-Location }

$hash = (Get-FileHash $bundle -Algorithm SHA256).Hash.ToLower()
Set-Content -Path $shaFile -Value "$hash  prp49_bundle.tar.gz" -Encoding ascii
$sizeMB = [math]::Round((Get-Item $bundle).Length / 1MB, 1)
Write-Host "DONE: $bundle ($sizeMB MB)"
Write-Host "SHA256: $hash"
Write-Host "下一步: scp prp49_bundle.tar.gz YOU@sylogin1.hpc.sjtu.edu.cn:~/"

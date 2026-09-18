param(
    [string]$OutputDir = (Join-Path $PSScriptRoot "native"),
    [string]$CudaArch = "all-major"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if (-not (Get-Command nvcc -ErrorAction SilentlyContinue)) {
    throw "nvcc was not found. Install a CUDA toolkit or build a CPU-only artifact."
}

if ($env:JAVA_HOME) {
    $Jdk = $env:JAVA_HOME
} else {
    $Javac = (Get-Command javac -ErrorAction Stop).Source
    $Jdk = Split-Path (Split-Path $Javac -Parent) -Parent
}

$JniInclude = Join-Path $Jdk "include"
$JniPlatformInclude = Join-Path $JniInclude "win32"
if (-not (Test-Path (Join-Path $JniPlatformInclude "jni_md.h"))) {
    throw "Could not find Windows JNI headers under $JniPlatformInclude"
}

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

$MinCudaCc = 0
if ($CudaArch -eq "all-major") {
    $Supported = & nvcc --list-gpu-arch |
        ForEach-Object { if ($_ -match '^compute_(\d+)$') { [int]$Matches[1] } } |
        Sort-Object
    if ($Supported) { $MinCudaCc = $Supported[0] }
} elseif ($CudaArch -match '^(?:sm|compute)_(\d+)$') {
    $MinCudaCc = [int]$Matches[1]
}

$Common = @(
    "-arch=$CudaArch",
    "-O3",
    "--shared",
    "-Xcompiler=/MD",
    "-I$JniInclude",
    "-I$JniPlatformInclude",
    "-DSTELARX_MIN_CUDA_CC=$MinCudaCc"
)

$Libraries = @(
    @{ Source = "stelarx_weight.cu";     Output = "stelarx_weight.dll" },
    @{ Source = "stelarx_dp.cu";         Output = "stelarx_dp.dll" },
    @{ Source = "stelarx_dist.cu";       Output = "stelarx_dist.dll" },
    @{ Source = "stelarx_similarity.cu"; Output = "stelarx_sim.dll" }
)

Write-Host "=== Building STELAR-X native GPU libraries ==="
Write-Host "  JDK         : $Jdk"
Write-Host "  CUDA arch   : $CudaArch"
Write-Host "  Minimum CC  : $MinCudaCc"
Write-Host "  Output      : $OutputDir"

foreach ($Library in $Libraries) {
    $Source = Join-Path $PSScriptRoot (Join-Path "src/native" $Library.Source)
    $Output = Join-Path $OutputDir $Library.Output
    Write-Host "  Building    : $Source -> $Output"
    & nvcc @Common "-o" $Output $Source
    if ($LASTEXITCODE -ne 0) {
        throw "nvcc failed while building $($Library.Source) (exit $LASTEXITCODE)"
    }
    Write-Host "  OK"
}

Write-Host "=== Native build complete ==="

param(
    [switch]$CpuOnly,
    [switch]$WithoutCuda,
    [switch]$WithCuda,
    [string]$CudaArch = "all-major",
    [string]$Version = "",
    [string]$OutputDir = (Join-Path $PSScriptRoot "dist"),
    [switch]$NoArchive,
    [switch]$Force
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$DisableCuda = $CpuOnly -or $WithoutCuda
if ($DisableCuda -and $WithCuda) {
    throw "Use either -WithoutCuda/-CpuOnly or -WithCuda, not both."
}

foreach ($Tool in @("javac", "jar", "jlink", "jpackage")) {
    if (-not (Get-Command $Tool -ErrorAction SilentlyContinue)) {
        throw "'$Tool' is required to build the artifact (JDK 21 or newer expected)."
    }
}

$JavacVersion = (& javac -version 2>&1).ToString()
if ($JavacVersion -notmatch '(\d+)') {
    throw "Could not determine javac version from '$JavacVersion'."
}
if ([int]$Matches[1] -lt 21) {
    throw "JDK 21 or newer is required; found $JavacVersion."
}

$VersionSource = Get-Content (Join-Path $PSScriptRoot "src/stelarx/Version.java") -Raw
if ($VersionSource -notmatch 'DEFAULT\s*=\s*"([^"]+)"') {
    throw "Could not determine the STELAR-X source version."
}
$SourceVersion = $Matches[1]
if (-not $Version) { $Version = $SourceVersion }
if ($Version -notmatch '^\d+(?:\.\d+){0,2}$') {
    throw "-Version must contain one to three numeric components (for example 1.2.0)."
}
$Arch = [System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture.ToString().ToLowerInvariant()
$PlatformArch = switch ($Arch) {
    "x64"   { "x86_64" }
    "arm64" { "arm64" }
    default { $Arch }
}

$NvccAvailable = [bool](Get-Command nvcc -ErrorAction SilentlyContinue)
$IncludeCuda = -not $DisableCuda -and $NvccAvailable
if ($WithCuda -and -not $NvccAvailable) {
    throw "-WithCuda was requested, but nvcc was not found."
}
if ($IncludeCuda -and $PlatformArch -ne "x86_64") {
    if ($WithCuda) {
        throw "The Windows CUDA artifact is currently supported only on x86-64."
    }
    $IncludeCuda = $false
}

$Capability = if ($IncludeCuda) { "cuda-with-cpu-fallback" } else { "cpu" }
# One public artifact per OS/CPU family. The CUDA-enabled image already carries
# the complete CPU implementation and falls back automatically.
$Artifact = "stelarx-$Version-windows-$PlatformArch"
$VersionDir = Join-Path $OutputDir $Version
$FinalImage = Join-Path $VersionDir $Artifact
$Archive = Join-Path $VersionDir "$Artifact.zip"
$ManifestPath = Join-Path $VersionDir "$Artifact.manifest.json"
if (-not $Force -and ((Test-Path $FinalImage) -or (Test-Path $Archive) -or
        (Test-Path $ManifestPath))) {
    throw "Artifact already exists for STELAR-X $Version on windows-$PlatformArch. " +
          "Use -Force to replace only this platform artifact, or choose another -Version."
}
$Work = Join-Path ([System.IO.Path]::GetTempPath()) ("stelarx-portable-" + [guid]::NewGuid())

try {
    New-Item -ItemType Directory -Force -Path $Work, $VersionDir | Out-Null
    $Build = Join-Path $Work "classes"
    $InputDir = Join-Path $Work "input"
    $Runtime = Join-Path $Work "runtime"
    $JpackageOut = Join-Path $Work "jpackage"
    New-Item -ItemType Directory -Force -Path $Build, $InputDir, $JpackageOut | Out-Null

    Write-Host "=== STELAR-X portable build ==="
    Write-Host "  Version      : $Version"
    Write-Host "  Platform     : windows-$PlatformArch"
    Write-Host "  CUDA bundle  : $IncludeCuda"
    Write-Host "  Artifact     : $Artifact"

    $Sources = Get-ChildItem (Join-Path $PSScriptRoot "src") -Recurse -Filter *.java |
        ForEach-Object { $_.FullName }
    & javac -d $Build -sourcepath (Join-Path $PSScriptRoot "src") @Sources
    if ($LASTEXITCODE -ne 0) { throw "javac failed (exit $LASTEXITCODE)" }

    $JarPath = Join-Path $InputDir "stelarx.jar"
    $JarManifest = Join-Path $Work "MANIFEST.MF"
    @(
        "Manifest-Version: 1.0",
        "Main-Class: stelarx.Main",
        "Implementation-Title: STELAR-X",
        "Implementation-Version: $Version",
        ""
    ) | Set-Content -Path $JarManifest -Encoding ASCII
    & jar --create --file $JarPath --manifest $JarManifest -C $Build .
    if ($LASTEXITCODE -ne 0) { throw "jar failed (exit $LASTEXITCODE)" }

    if ($IncludeCuda) {
        & (Join-Path $PSScriptRoot "build_native.ps1") -OutputDir $InputDir -CudaArch $CudaArch
    }

    & jlink --add-modules java.base --strip-debug --no-header-files --no-man-pages `
        --compress=zip-6 --output $Runtime
    if ($LASTEXITCODE -ne 0) { throw "jlink failed (exit $LASTEXITCODE)" }

    & jpackage --type app-image --name stelarx --app-version $Version `
        --input $InputDir --main-jar stelarx.jar --main-class stelarx.Main `
        --runtime-image $Runtime --dest $JpackageOut `
        --java-options '-Djava.library.path=$APPDIR' `
        --java-options '-Dfile.encoding=UTF-8' `
        --java-options '-XX:InitialRAMPercentage=2.0' `
        --java-options '-XX:MaxRAMPercentage=85.0' `
        --java-options '-XX:ErrorFile=crash_logs/stelarx-hotspot-crash-%p.log'
    if ($LASTEXITCODE -ne 0) { throw "jpackage failed (exit $LASTEXITCODE)" }

    $Image = Join-Path $JpackageOut "stelarx"
    $ExampleDir = Join-Path $Image "example"
    New-Item -ItemType Directory -Force -Path $ExampleDir | Out-Null
    Copy-Item (Join-Path $PSScriptRoot "all_gt_bs_rooted_37.tre") `
        (Join-Path $ExampleDir "all_gt_37.tre")
    Copy-Item (Join-Path $PSScriptRoot "true_37.tre") `
        (Join-Path $ExampleDir "true_37.tre")

    $Readme = @"
STELAR-X $Version - self-contained windows-$PlatformArch build

Run in PowerShell or Command Prompt:
  .\stelarx.exe --help
  .\stelarx.exe --diagnose
  .\stelarx.exe -i C:\path\to\rooted_gene_trees.tre -o C:\path\to\output_species_tree.tre

Ready-made 37-taxon example (run from this directory):
  .\stelarx.exe -i example\all_gt_37.tre -o example\predicted_st_37.tre --search-space S1 -vv

The example directory initially contains:
  all_gt_37.tre   input gene trees
  true_37.tre     reference/true species tree

After the command finishes, example\predicted_st_37.tre contains the inferred
species tree. The reference tree is provided for comparison and is not used by
STELAR-X during inference.

No Java installation is required. This artifact includes CUDA acceleration: $IncludeCuda.
CUDA still requires a compatible NVIDIA GPU and installed NVIDIA driver. If CUDA cannot
be used, STELAR-X automatically explains why and falls back to CPU. Use --gpu-strict only
when falling back would be undesirable.

Unexpected failure reports are stored in crash_logs under the directory from
which you launch STELAR-X. Set STELARX_CRASH_DIR to override Java report storage.
"@
    Set-Content -Path (Join-Path $Image "README.txt") -Value $Readme -Encoding UTF8

    $BuildInfo = @(
        "version=$Version",
        "platform=windows-$PlatformArch",
        "capability=$Capability",
        "built_utc=$([DateTime]::UtcNow.ToString('yyyy-MM-ddTHH:mm:ssZ'))",
        "java=$((& java -version 2>&1 | Select-Object -First 1).ToString())"
    )
    if ($IncludeCuda) {
        $BuildInfo += "nvcc=$((& nvcc --version | Select-Object -Last 1).ToString())"
        $BuildInfo += "cuda_arch=$CudaArch"
    }
    if (Get-Command git -ErrorAction SilentlyContinue) {
        $Commit = (& git -C $PSScriptRoot rev-parse HEAD 2>$null)
        if (-not $Commit) { $Commit = "unknown" }
        $BuildInfo += "git_commit=$Commit"
    }
    Set-Content -Path (Join-Path $Image "BUILD-INFO.txt") -Value $BuildInfo -Encoding UTF8

    $Launcher = Join-Path $Image "stelarx.exe"
    $PreviousNoColor = $env:NO_COLOR
    try {
        $env:NO_COLOR = "1"
        $VersionOutput = (& $Launcher --version | Out-String).Trim()
        $VersionExitCode = $LASTEXITCODE
    }
    finally {
        if ($null -eq $PreviousNoColor) {
            Remove-Item Env:NO_COLOR -ErrorAction SilentlyContinue
        }
        else {
            $env:NO_COLOR = $PreviousNoColor
        }
    }
    if ($VersionExitCode -ne 0 `
            -or -not $VersionOutput.Contains("STELAR-X  v$Version") `
            -or -not $VersionOutput.Contains("Welcome to STELAR-X version $Version!")) {
        throw "Packaged --version smoke test failed: '$VersionOutput'."
    }
    & $Launcher --cpu --diagnose | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "Packaged --diagnose smoke test failed." }
    $SmokeTree = Join-Path $Work "smoke-species-tree.tre"
    & $Launcher --cpu --search-space S2 -q `
        -i (Join-Path $ExampleDir "all_gt_37.tre") -o $SmokeTree `
        --log-file (Join-Path $Work "smoke-run.log")
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path $SmokeTree) -or
            (Get-Item $SmokeTree).Length -eq 0) {
        throw "Packaged end-to-end inference smoke test failed."
    }
    $SmokeLog = Join-Path $Work "smoke-run.log"
    if (-not (Test-Path $SmokeLog) -or
            -not (Select-String -Path $SmokeLog -SimpleMatch "Run Summary" -Quiet)) {
        throw "Packaged --log-file smoke test did not capture the complete run."
    }

    if (Test-Path $FinalImage) { Remove-Item -Recurse -Force $FinalImage }
    Move-Item $Image $FinalImage

    $ArchiveName = ""
    $ArchiveHash = ""
    if (-not $NoArchive) {
        if (Test-Path $Archive) { Remove-Item -Force $Archive }
        Compress-Archive -Path $FinalImage -DestinationPath $Archive -CompressionLevel Optimal
        $ArchiveHash = (Get-FileHash -Algorithm SHA256 $Archive).Hash.ToLowerInvariant()
        $ArchiveName = [IO.Path]::GetFileName($Archive)
        Set-Content -Path "$Archive.sha256" -Value "$ArchiveHash  $ArchiveName" -Encoding ASCII
        Write-Host "  Archive      : $Archive"
        Write-Host "  SHA-256      : $Archive.sha256"
    }

    [ordered]@{
        version = $Version
        platform = "windows-$PlatformArch"
        capability = $Capability
        minimum_glibc = ""
        archive = $ArchiveName
        sha256 = $ArchiveHash
    } | ConvertTo-Json | Set-Content -Path $ManifestPath -Encoding UTF8

    Write-Host "Portable application ready: $FinalImage\stelarx.exe"
    Write-Host "Release manifest: $ManifestPath"
}
finally {
    if (Test-Path $Work) { Remove-Item -Recurse -Force $Work }
}

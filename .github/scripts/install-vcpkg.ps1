param(
    [Parameter(Mandatory = $true)][string]$Triplet,
    [Parameter(Mandatory = $true)][string]$VcpkgRoot
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "../..")).Path
$cacheDir = Join-Path $repoRoot ".vcpkg-cache"
$vcpkgTag = "2026.07.29"

New-Item -ItemType Directory -Force -Path $cacheDir | Out-Null
if (-not (Test-Path (Join-Path $VcpkgRoot ".git"))) {
    git clone --branch $vcpkgTag --depth 1 https://github.com/microsoft/vcpkg.git $VcpkgRoot
    if ($LASTEXITCODE -ne 0) { throw "vcpkg clone failed" }
}

& (Join-Path $VcpkgRoot "bootstrap-vcpkg.bat") -disableMetrics
if ($LASTEXITCODE -ne 0) { throw "vcpkg bootstrap failed" }

$env:VCPKG_BINARY_SOURCES = "clear;files,$cacheDir,readwrite"
& (Join-Path $VcpkgRoot "vcpkg.exe") install "flint:$Triplet" "--overlay-triplets=$repoRoot/.github/triplets"
if ($LASTEXITCODE -ne 0) { throw "vcpkg dependency installation failed" }

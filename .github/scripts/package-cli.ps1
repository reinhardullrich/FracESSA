param(
    [Parameter(Mandatory = $true)][string]$Version,
    [Parameter(Mandatory = $true)][string]$Platform,
    [Parameter(Mandatory = $true)][string]$Triplet,
    [Parameter(Mandatory = $true)][string]$VcpkgRoot
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "../..")).Path
$buildDir = Join-Path $repoRoot "cpp/build-release-$Platform"
$output = Join-Path $repoRoot "dist/fracessa-$Version-$Platform.exe"
$tripletsDir = Join-Path $repoRoot ".github/triplets"

& (Join-Path $PSScriptRoot "install-vcpkg.ps1") -Triplet $Triplet -VcpkgRoot $VcpkgRoot

cmake -B $buildDir -S (Join-Path $repoRoot "cpp") `
    -DCMAKE_BUILD_TYPE=Release `
    "-DCMAKE_TOOLCHAIN_FILE=$VcpkgRoot/scripts/buildsystems/vcpkg.cmake" `
    "-DVCPKG_TARGET_TRIPLET=$Triplet" `
    "-DVCPKG_OVERLAY_TRIPLETS=$tripletsDir" `
    -DBUILD_TESTING=ON `
    -DFRACESSA_NATIVE_ARCH=OFF `
    -DFRACESSA_BUILD_CLI=ON `
    -DFRACESSA_BUILD_PYTHON=OFF
if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }

cmake --build $buildDir --config Release -j 4
if ($LASTEXITCODE -ne 0) { throw "CMake build failed" }

ctest --test-dir $buildDir --build-config Release --output-on-failure --parallel 4
if ($LASTEXITCODE -ne 0) { throw "CTest failed" }

$binary = Join-Path $buildDir "Release/fracessa.exe"
$actualVersion = (& $binary --version).Trim()
if ($actualVersion -ne $Version) {
    throw "Release version $Version does not match binary version $actualVersion"
}

New-Item -ItemType Directory -Force -Path (Split-Path $output) | Out-Null
Copy-Item $binary $output -Force

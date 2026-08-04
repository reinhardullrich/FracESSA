param(
    [Parameter(Mandatory = $true)][string]$Version,
    [Parameter(Mandatory = $true)][string]$Platform,
    [Parameter(Mandatory = $true)][string]$Triplet,
    [Parameter(Mandatory = $true)][string]$VcpkgRoot
)

$ErrorActionPreference = "Stop"
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "../..")).Path
$buildDir = Join-Path $repoRoot "cpp/build-release-$Platform"
$packageName = "fracessa-$Version-$Platform"
$packageDir = Join-Path $repoRoot "dist/$packageName"
$tripletsDir = Join-Path $repoRoot ".github/triplets"

& (Join-Path $PSScriptRoot "install-vcpkg.ps1") -Triplet $Triplet -VcpkgRoot $VcpkgRoot

cmake -B $buildDir -S (Join-Path $repoRoot "cpp") `
    -DCMAKE_BUILD_TYPE=Release `
    "-DCMAKE_TOOLCHAIN_FILE=$VcpkgRoot/scripts/buildsystems/vcpkg.cmake" `
    "-DVCPKG_TARGET_TRIPLET=$Triplet" `
    "-DVCPKG_OVERLAY_TRIPLETS=$tripletsDir" `
    -DBUILD_TESTING=OFF `
    -DFRACESSA_NATIVE_ARCH=OFF `
    -DFRACESSA_BUILD_CLI=ON `
    -DFRACESSA_BUILD_PYTHON=OFF
if ($LASTEXITCODE -ne 0) { throw "CMake configure failed" }

cmake --build $buildDir --config Release --target fracessa -j 4
if ($LASTEXITCODE -ne 0) { throw "CMake build failed" }

$binary = Join-Path $buildDir "Release/fracessa.exe"
$actualVersion = (& $binary --version).Trim()
if ($actualVersion -ne $Version) {
    throw "Tag version $Version does not match binary version $actualVersion"
}

New-Item -ItemType Directory -Force -Path $packageDir | Out-Null
Copy-Item $binary (Join-Path $packageDir "fracessa.exe")
Copy-Item (Join-Path $repoRoot "README.md") (Join-Path $packageDir "README.md")
Copy-Item (Join-Path $repoRoot "LICENSE") (Join-Path $packageDir "LICENSE")
Copy-Item (Join-Path $repoRoot "THIRD_PARTY_NOTICES.md") (Join-Path $packageDir "THIRD_PARTY_NOTICES.md")
Compress-Archive -Path $packageDir -DestinationPath (Join-Path $repoRoot "dist/$packageName.zip") -Force

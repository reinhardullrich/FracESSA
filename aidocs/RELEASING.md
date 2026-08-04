# Releasing FracESSA

Ordinary pushes and pull requests run no GitHub Actions. A release starts only when a version tag matching `v*` is pushed.

1. Set the calendar version in `cpp/CMakeLists.txt` using `YYYY.M.D.N`, where `N` starts at 1 and increments for additional releases on the same day.
2. Commit and push the release-ready source.
3. Create and push the matching tag, for example `v2026.8.4.1`.

The tag and CMake version must match exactly. The release workflow then builds and tests:

- standalone CLI archives for Linux x86-64, Linux ARM64, macOS Intel, macOS Apple Silicon, and Windows x86-64;
- `pyfracessa` wheels for CPython 3.11, 3.12, 3.13, and 3.14 on those five platforms;
- one Python source distribution.

The wheels and source distribution are published to PyPI through the `pypi` GitHub environment and PyPI trusted publishing. No API token is stored. The five CLI archives are attached to the GitHub release. Linux artifacts target glibc 2.28 or newer, and release builds disable `-march=native`.

FracESSA and the complete statically linked release are distributed under GPL-3.0-or-later. Every CLI archive and Python wheel
contains the project license and THIRD_PARTY_NOTICES.md, which identifies the exact corresponding source and vcpkg patches.

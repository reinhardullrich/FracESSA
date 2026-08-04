# Releasing FracESSA

Ordinary pushes and pull requests run no GitHub Actions. A release starts only when the `Release` workflow is manually run from `main`.

1. Set the calendar version in `cpp/CMakeLists.txt` using `YYYY.M.D.N`, where `N` starts at 1 and increments for additional releases on the same day.
2. Commit and push the release-ready source.
3. In GitHub Actions, open `Release`, choose `Run workflow`, select `main`, and start it.

The workflow refuses to run from another branch or to reuse an existing version tag. It builds and tests:

- standalone CLI archives for Linux x86-64, Linux ARM64, macOS Intel, macOS Apple Silicon, and Windows x86-64;
- `pyfracessa` wheels for CPython 3.11, 3.12, 3.13, and 3.14 on those five platforms;
- one Python source distribution.

Only after every build succeeds does the workflow create the matching tag and GitHub release. It attaches the five CLI archives, then publishes the wheels and source distribution to PyPI through the `pypi` GitHub environment and PyPI trusted publishing. A failed build therefore creates neither a tag nor a release. No API token is stored. Linux artifacts target glibc 2.28 or newer, and release builds disable `-march=native`.

Release jobs run from `main`, so their vcpkg binary caches are shared by later releases. GMP, MPFR, and FLINT rebuild only when their versions or release triplets change.

FracESSA and the complete statically linked release are distributed under GPL-3.0-or-later. Every CLI archive and Python wheel
contains the project license and THIRD_PARTY_NOTICES.md, which identifies the exact corresponding source and vcpkg patches.

# Third-party notices and corresponding source

FracESSA is licensed under GPL-3.0-or-later. Prebuilt releases statically link third-party software.

For a FracESSA binary reporting version V, its complete project source, build scripts, and release configuration are in the
matching vV tag at <https://github.com/reinhardullrich/fracessa/tags>. GitHub also provides source archives beside every release.

The exact dependency sources used by the release are available without charge here:

| Component | Version | License | Exact source |
|---|---:|---|---|
| FLINT | 3.6.0 | LGPL-3.0-or-later | <https://github.com/flintlib/flint/archive/refs/tags/v3.6.0.tar.gz> |
| GMP | 6.3.0 | LGPL-3.0-only OR GPL-2.0-only | <https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz> |
| MPFR | 4.2.2 | LGPL-3.0-or-later | <https://www.mpfr.org/mpfr-4.2.2/mpfr-4.2.2.tar.xz> |
| spdlog | 1.13.0 | MIT | <https://github.com/gabime/spdlog/archive/refs/tags/v1.13.0.tar.gz> |
| fmt, bundled by spdlog | 9.1.0 | MIT | Included in the spdlog source above |
| argparse | 2.9 | MIT | <https://github.com/p-ranav/argparse/archive/refs/tags/v2.9.tar.gz> |
| pybind11 | 3.0.4 | BSD-3-Clause | <https://github.com/pybind/pybind11/archive/refs/tags/v3.0.4.tar.gz> |
| PThreads4W, Windows only | 3.0.0 | Apache-2.0 | <https://sourceforge.net/projects/pthreads4w/files/pthreads4w-code-v3.0.0.zip/download> |
| gettimeofday compatibility code | 2017-10-14 | Permissive PostgreSQL-derived license | See the pinned vcpkg port below |

The release uses vcpkg commit 9e593bb18ea69cc5095e012465dcd675a822ed0d. Its exact recipes and patches are at:

- <https://github.com/microsoft/vcpkg/tree/9e593bb18ea69cc5095e012465dcd675a822ed0d/ports/flint>
- <https://github.com/microsoft/vcpkg/tree/9e593bb18ea69cc5095e012465dcd675a822ed0d/ports/gmp>
- <https://github.com/microsoft/vcpkg/tree/9e593bb18ea69cc5095e012465dcd675a822ed0d/ports/mpfr>
- <https://github.com/microsoft/vcpkg/tree/9e593bb18ea69cc5095e012465dcd675a822ed0d/ports/pthreads>
- <https://github.com/microsoft/vcpkg/tree/9e593bb18ea69cc5095e012465dcd675a822ed0d/ports/gettimeofday>

Each linked source archive contains its copyright and license text. The relevant online license notices are:

- GNU GPLv3: <https://www.gnu.org/licenses/gpl-3.0.html>
- GNU LGPLv3: <https://www.gnu.org/licenses/lgpl-3.0.html>
- Apache License 2.0: <https://www.apache.org/licenses/LICENSE-2.0>
- spdlog and bundled fmt: <https://github.com/gabime/spdlog/blob/v1.13.0/LICENSE>
- argparse: <https://github.com/p-ranav/argparse/blob/v2.9/LICENSE>
- pybind11: <https://github.com/pybind/pybind11/blob/v3.0.4/LICENSE>

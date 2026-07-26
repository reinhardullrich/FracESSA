#!/usr/bin/env python3
"""Make archived tags use the common ARM64 math libraries without changing code."""

from pathlib import Path
import re
import sys


path = Path(sys.argv[1])
text = path.read_text()
if "ExternalProject_Add(flint_build" not in text:
    raise SystemExit(0)

replacement = """else()
    # Experiment-only build adaptation: use one common native FLINT installation.
    find_path(FLINT_INCLUDE_DIR NAMES flint/flint.h)
    find_library(FLINT_LIB NAMES flint)
    if(FLINT_LIB AND FLINT_INCLUDE_DIR)
        set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})
        set(FLINT_LIBRARIES ${FLINT_LIB})
    else()
        message(FATAL_ERROR "Native FLINT not found for historical benchmark")
    endif()
endif()"""

pattern = re.compile(
    r"else\(\)\n\s*# LINUX / MACOS:.*?ExternalProject_Add\(flint_build.*?\nendif\(\)",
    re.DOTALL,
)
text, count = pattern.subn(replacement, text, count=1)
if count != 1:
    raise RuntimeError(f"could not replace FLINT ExternalProject block in {path}")

text = re.sub(r"^\s*add_dependencies\(fracessa_lib flint_build\)\s*$", "", text, flags=re.MULTILINE)
path.write_text(text)


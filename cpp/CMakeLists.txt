cmake_minimum_required(VERSION 3.18)

project(fracessa
    VERSION 1.0.0
    DESCRIPTION "FRACESSA - Fractional ESS Analyzer"
    LANGUAGES CXX C
)

# ==============================================================================
# 1. Compiler Settings & Optimization
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
option(FRACESSA_NATIVE_ARCH "Enable local CPU-specific optimizations with -march=native" ON)

if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    add_compile_options(/O2 /Oi /Ot /GL /fp:fast /DNDEBUG)
    add_link_options(/LTCG)
    # Suppress noisy warnings from GMP/FLINT headers on Windows
    add_compile_options(/wd4244 /wd4146 /wd4101)
    # Ensure the linker finds all static dependencies correctly
    add_link_options(/NODEFAULTLIB:MSVCRT /NODEFAULTLIB:MSVCRTD)
else()
    add_compile_options(-O3 -funroll-loops -DNDEBUG)
    if(FRACESSA_NATIVE_ARCH)
        add_compile_options(-march=native)
    endif()

    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED)
    if(IPO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()

# ==============================================================================
# 2. Dependencies (spdlog, argparse)
# ==============================================================================
include(FetchContent)

FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.13.0
)
FetchContent_Declare(argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v2.9
)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v3.0.4
)

# Use the same `python3` selected by project scripts unless explicitly overridden.
if(NOT Python3_EXECUTABLE)
    find_program(Python3_EXECUTABLE NAMES python3 python REQUIRED)
endif()
find_package(Python3 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
set(PYBIND11_FINDPYTHON ON CACHE BOOL "Use CMake FindPython for pybind11" FORCE)

# MSVC-specific settings for FetchContent
if(MSVC)
    set(gtest_force_shared_crt OFF CACHE BOOL "" FORCE)
endif()

FetchContent_MakeAvailable(spdlog argparse googletest pybind11)

enable_testing()

# ==============================================================================
# 3. Math Dependencies (GMP, MPFR, FLINT)
# ==============================================================================

# --- GMP & MPFR (Generic search for all platforms) ---
if(UNIX AND NOT APPLE)
    # Shared libs are required for linking the pybind module.
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
elseif(APPLE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib" ".a")
endif()

# Search for GMP (Common names: gmp, libgmp, or mpir on Windows)
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIB NAMES gmp libgmp mpir libmpir)

# Search for MPFR
find_path(MPFR_INCLUDE_DIR NAMES mpfr.h)
find_library(MPFR_LIB NAMES mpfr libmpfr)

if(NOT (GMP_LIB AND GMP_INCLUDE_DIR))
    message(FATAL_ERROR "GMP not found.")
endif()
if(NOT (MPFR_LIB AND MPFR_INCLUDE_DIR))
    message(FATAL_ERROR "MPFR not found.")
endif()

# --- FLINT ---
set(FLINT_INCLUDE_DIRS "")
set(FLINT_LIBRARIES "")

if(WIN32)
    # WINDOWS (vcpkg): Search for pre-installed static FLINT
    find_path(FLINT_INCLUDE_DIR NAMES flint/flint.h)
    find_library(FLINT_LIB NAMES flint libflint)

    if(FLINT_LIB AND FLINT_INCLUDE_DIR)
        set(FLINT_LIBRARIES ${FLINT_LIB})
        set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})

        # Static FLINT on Windows often needs pthreads
        find_library(PTHREADS_LIB NAMES pthreadVC3 pthreadVC2 pthread pthreads)
        if(PTHREADS_LIB)
             list(APPEND FLINT_LIBRARIES ${PTHREADS_LIB})
        endif()
    else()
        message(FATAL_ERROR "FLINT static library not found.")
    endif()
else()
    # LINUX / MACOS: Use system-installed FLINT
    find_path(FLINT_INCLUDE_DIR NAMES flint/flint.h)
    find_library(FLINT_LIB NAMES flint)

    if(FLINT_LIB AND FLINT_INCLUDE_DIR)
        set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})
        set(FLINT_LIBRARIES ${FLINT_LIB})
        message(STATUS "Found system FLINT: ${FLINT_LIB}")
    else()
        message(FATAL_ERROR "FLINT library not found. Install it with 'sudo apt install libflint-dev' or 'brew install flint'")
    endif()
endif()

# ==============================================================================
# 4. Build Targets
# ==============================================================================

add_library(fracessa_find_candidate_verified OBJECT
    src/find_candidate_verified.cpp
)
set_target_properties(fracessa_find_candidate_verified PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION FALSE
    POSITION_INDEPENDENT_CODE ON
)
target_include_directories(fracessa_find_candidate_verified PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fracessa
    ${FLINT_INCLUDE_DIRS}
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
)

if(MSVC)
    target_compile_options(fracessa_find_candidate_verified PRIVATE /fp:strict /GL-)
    set(FRACESSA_FIND_CANDIDATE_VERIFIED_SUPPORTED 1)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(fracessa_find_candidate_verified PRIVATE -fno-fast-math -ffp-contract=off)
    set(FRACESSA_FIND_CANDIDATE_VERIFIED_SUPPORTED 1)
else()
    set(FRACESSA_FIND_CANDIDATE_VERIFIED_SUPPORTED 0)
endif()
target_compile_definitions(fracessa_find_candidate_verified PRIVATE
    FRACESSA_FIND_CANDIDATE_VERIFIED_SUPPORTED=${FRACESSA_FIND_CANDIDATE_VERIFIED_SUPPORTED}
)

add_library(fracessa_lib
    src/fracessa.cpp
    src/find_candidate_unsafe.cpp
    src/find_candidate_exact.cpp
    src/checkstab.cpp
    $<TARGET_OBJECTS:fracessa_find_candidate_verified>
)
set_target_properties(fracessa_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Set common include directories
target_include_directories(fracessa_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/fracessa>
    ${FLINT_INCLUDE_DIRS}
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
)

# Link all math libraries. Order matters for static linking: FLINT -> MPFR -> GMP
target_link_libraries(fracessa_lib PUBLIC
    spdlog::spdlog
    ${FLINT_LIBRARIES}
    ${MPFR_LIB}
    ${GMP_LIB}
)


add_executable(fracessa src/main.cpp)

target_link_libraries(fracessa PRIVATE
    fracessa_lib
    argparse::argparse
)

pybind11_add_module(fracessa_core MODULE src/pybind_module.cpp)
target_link_libraries(fracessa_core PRIVATE fracessa_lib)

install(TARGETS fracessa RUNTIME DESTINATION bin)

add_subdirectory(tests)

cmake_minimum_required(VERSION 3.14)

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

if(MSVC)
    add_compile_options(/O2 /Oi /Ot /GL /fp:fast /DNDEBUG)
    add_link_options(/LTCG)
    # Suppress noisy warnings from GMP/FLINT headers on Windows
    add_compile_options(/wd4244 /wd4146 /wd4101)
    # Ensure the linker finds all static dependencies correctly
    add_link_options(/NODEFAULTLIB:library)
else()
    add_compile_options(-O3 -march=native -flto -funroll-loops -DNDEBUG)

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
FetchContent_MakeAvailable(spdlog argparse googletest)

enable_testing()

# ==============================================================================
# 3. Math Dependencies (GMP, MPFR, FLINT)
# ==============================================================================

# --- GMP & MPFR (Generic search for all platforms) ---
if(NOT WIN32)
    # Force static suffixes for Linux/Mac
    set(_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# Search for GMP (Common names: gmp, libgmp, or mpir on Windows)
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIB NAMES gmp libgmp mpir libmpir)

# Search for MPFR
find_path(MPFR_INCLUDE_DIR NAMES mpfr.h)
find_library(MPFR_LIB NAMES mpfr libmpfr)

if(NOT WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

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

add_library(fracessa_lib
    src/fracessa.cpp
    src/findeq.cpp
    src/checkstab.cpp
)

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

if(NOT WIN32)
    find_package(Threads REQUIRED)
    target_link_libraries(fracessa_lib PUBLIC Threads::Threads)
endif()

add_executable(fracessa src/main.cpp)

target_link_libraries(fracessa PRIVATE
    fracessa_lib
    argparse::argparse
)

install(TARGETS fracessa RUNTIME DESTINATION bin)

add_subdirectory(tests)
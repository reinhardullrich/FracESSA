cmake_minimum_required(VERSION 3.14)

project(fracessa
    VERSION 1.0.0
    DESCRIPTION "FRACESSA - Fractional ESS Analyzer"
    LANGUAGES CXX C
)

# ==============================================================================
# 1. Compiler Settings
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    # Windows: VCPKG static triplet handles /MT (static runtime) automatically
    add_compile_options(/O2 /Oi /Ot /GL /fp:fast /DNDEBUG)
    add_link_options(/LTCG)
    # Suppress warnings about unreferenced variables and conversion loss (common in GMP/FLINT headers)
    add_compile_options(/wd4244 /wd4146 /wd4101) 
else()
    # Linux/Mac: Optimization flags
    add_compile_options(-O3 -march=native -flto -funroll-loops -DNDEBUG)
    
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        add_link_options(-static-libgcc -static-libstdc++)
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
FetchContent_MakeAvailable(spdlog argparse)

# ==============================================================================
# 3. Math Dependencies (GMP, MPFR, FLINT)
# ==============================================================================

if(WIN32)
    # --- WINDOWS (VCPKG) ---
    # Use find_package which works reliably with vcpkg
    # Vcpkg provides "Unofficial" configs for GMP/MPFR usually
    find_package(gmp REQUIRED)
    find_package(mpfr REQUIRED)
    find_package(flint REQUIRED)
    
    # Map vcpkg targets to variables for consistency
    # Note: different vcpkg versions might expose different target names.
    # Common ones: gmp::gmp, mpfr::mpfr, flint::flint
    set(GMP_LIBRARIES gmp::gmp)
    set(MPFR_LIBRARIES mpfr::mpfr)
    set(FLINT_LIBRARIES flint::flint)
    
    # FLINT often needs pthread on Windows if built with thread support
    find_package(PThreads)
    if(PThreads_FOUND)
        list(APPEND FLINT_LIBRARIES PThreads::PThreads)
    endif()
    
    message(STATUS "Found FLINT (Windows vcpkg)")

else()
    # --- LINUX / MACOS (Manual / Source) ---
    
    # Force Static Linking on Linux/Mac
    set(_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    
    find_library(GMP_LIB NAMES gmp libgmp mpir)
    find_path(GMP_INCLUDE_DIR NAMES gmp.h)

    find_library(MPFR_LIB NAMES mpfr libmpfr)
    find_path(MPFR_INCLUDE_DIR NAMES mpfr.h)

    # Restore suffixes
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

    if(NOT (GMP_LIB AND GMP_INCLUDE_DIR))
        message(FATAL_ERROR "GMP static lib not found.")
    endif()

    if(NOT (MPFR_LIB AND MPFR_INCLUDE_DIR))
        message(FATAL_ERROR "MPFR static lib not found.")
    endif()
    
    # Build FLINT from source
    include(ExternalProject)
    get_filename_component(GMP_LIB_DIR ${GMP_LIB} DIRECTORY)
    get_filename_component(MPFR_LIB_DIR ${MPFR_LIB} DIRECTORY)

    ExternalProject_Add(flint_build
        URL https://github.com/flintlib/flint/releases/download/v3.0.1/flint-3.0.1.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/third_party/flint
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ./configure 
            --prefix=<INSTALL_DIR> 
            --with-gmp=${GMP_LIB_DIR}/.. 
            --with-mpfr=${MPFR_LIB_DIR}/..
            --disable-shared 
            --enable-static
        BUILD_COMMAND make -j4
        INSTALL_COMMAND make install
    )

    ExternalProject_Get_Property(flint_build INSTALL_DIR)
    set(FLINT_INCLUDE_DIRS ${INSTALL_DIR}/include)
    set(FLINT_LIBRARIES ${INSTALL_DIR}/lib/libflint.a)
    
    # Set variables for consistent usage below
    set(GMP_LIBRARIES ${GMP_LIB})
    set(MPFR_LIBRARIES ${MPFR_LIB})
    set(GMP_INCLUDE_DIR ${GMP_INCLUDE_DIR}) # Already set but consistent naming
    set(MPFR_INCLUDE_DIR ${MPFR_INCLUDE_DIR})
endif()

# ==============================================================================
# 4. Build Targets
# ==============================================================================

add_library(fracessa_lib
    src/fracessa.cpp
    src/findeq.cpp
    src/checkstab.cpp
)

# Standard Include Paths
target_include_directories(fracessa_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/fracessa>
)

# Math Include Paths (Only needed explicitly for non-vcpkg/source builds usually, but safer to add)
if(NOT WIN32)
    target_include_directories(fracessa_lib PUBLIC
        ${FLINT_INCLUDE_DIRS}
        ${GMP_INCLUDE_DIR}
        ${MPFR_INCLUDE_DIR}
    )
endif()

# Linking
target_link_libraries(fracessa_lib PUBLIC
    spdlog::spdlog
    ${FLINT_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${GMP_LIBRARIES}
)

if(NOT WIN32)
    add_dependencies(fracessa_lib flint_build)
    find_package(Threads REQUIRED)
    target_link_libraries(fracessa_lib PUBLIC Threads::Threads)
endif()

# Executable
add_executable(fracessa src/main.cpp)

target_link_libraries(fracessa PRIVATE
    fracessa_lib
    argparse::argparse
)

install(TARGETS fracessa RUNTIME DESTINATION bin)
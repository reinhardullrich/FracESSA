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
    add_compile_options(/O2 /Oi /Ot /GL /fp:fast /DNDEBUG)
    add_link_options(/LTCG)
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
    GIT_PROGRESS TRUE     # <--- Shows progress
)
FetchContent_Declare(argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v2.9
    GIT_PROGRESS TRUE     # <--- Shows progress
)
FetchContent_MakeAvailable(spdlog argparse)

# ==============================================================================
# 3. Math Dependencies (GMP, MPFR, FLINT)
# ==============================================================================

# --- GMP ---
find_library(GMP_LIB NAMES gmp libgmp mpir)
find_path(GMP_INCLUDE_DIR NAMES gmp.h)

if(NOT (GMP_LIB AND GMP_INCLUDE_DIR))
    message(FATAL_ERROR "GMP not found. \nLinux: sudo apt install libgmp-dev\nMac: brew install gmp\nWindows: vcpkg install gmp")
endif()

# --- MPFR ---
find_library(MPFR_LIB NAMES mpfr libmpfr)
find_path(MPFR_INCLUDE_DIR NAMES mpfr.h)

if(NOT (MPFR_LIB AND MPFR_INCLUDE_DIR))
    message(FATAL_ERROR "MPFR not found. \nLinux: sudo apt install libmpfr-dev\nMac: brew install mpfr\nWindows: vcpkg install mpfr")
endif()

# --- FLINT ---
set(FLINT_INCLUDE_DIRS "")
set(FLINT_LIBRARIES "")

if(WIN32)
    # WINDOWS: Use vcpkg provided libraries
    find_library(FLINT_LIB NAMES flint libflint)
    find_path(FLINT_INCLUDE_DIR NAMES flint/flint.h)
    
    if(FLINT_LIB AND FLINT_INCLUDE_DIR)
        set(FLINT_LIBRARIES ${FLINT_LIB})
        set(FLINT_INCLUDE_DIRS ${FLINT_INCLUDE_DIR})
        message(STATUS "Found FLINT (Windows): ${FLINT_LIB}")
    else()
        message(FATAL_ERROR "FLINT not found. Ensure 'vcpkg install flint:x64-windows' ran.")
    endif()
else()
    # LINUX / MACOS: Build FLINT 3.0.1 from Source (Tarball)
    # Using tarball ensures 'configure' script is present (fixes 'configure not found' error)
    include(ExternalProject)
    
    # Get directories to pass to FLINT configure
    get_filename_component(GMP_LIB_DIR ${GMP_LIB} DIRECTORY)
    get_filename_component(MPFR_LIB_DIR ${MPFR_LIB} DIRECTORY)

    ExternalProject_Add(flint_build
        URL https://github.com/flintlib/flint/releases/download/v3.0.1/flint-3.0.1.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/third_party/flint
        BUILD_IN_SOURCE 1
        # Configure FLINT to use the system's GMP/MPFR
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
endif()

# ==============================================================================
# 4. Build Targets
# ==============================================================================

# --- Library ---
add_library(fracessa_lib
    src/fracessa.cpp
    src/findeq.cpp
    src/checkstab.cpp
)

target_include_directories(fracessa_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/fracessa>
    ${FLINT_INCLUDE_DIRS}
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
)

target_link_libraries(fracessa_lib PUBLIC
    spdlog::spdlog
    ${FLINT_LIBRARIES}
    ${MPFR_LIB}
    ${GMP_LIB}
)

# On Linux/Mac, ensure FLINT is built before compiling the library
if(NOT WIN32)
    add_dependencies(fracessa_lib flint_build)
    find_package(Threads REQUIRED)
    target_link_libraries(fracessa_lib PUBLIC Threads::Threads)
endif()

# --- Executable ---
add_executable(fracessa src/main.cpp)

target_link_libraries(fracessa PRIVATE
    fracessa_lib
    argparse::argparse
)

# ==============================================================================
# 5. Installation
# ==============================================================================
install(TARGETS fracessa RUNTIME DESTINATION bin)
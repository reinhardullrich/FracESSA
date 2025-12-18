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
else()
    # Linux/Mac: Optimization flags
    add_compile_options(-O3 -march=native -flto -funroll-loops -DNDEBUG)
    
    # Linux: Link standard C++ libraries statically for maximum portability
    # This ensures the binary runs on older Linux distros without "GLIBCXX not found" errors
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

# --- GMP & MPFR (Linux/Mac Only) ---
# Windows handles this via vcpkg
if(NOT WIN32)
    # FORCE STATIC linking on Linux/Mac by looking for .a files only
    set(_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    
    find_library(GMP_LIB NAMES gmp libgmp mpir)
    find_path(GMP_INCLUDE_DIR NAMES gmp.h)

    find_library(MPFR_LIB NAMES mpfr libmpfr)
    find_path(MPFR_INCLUDE_DIR NAMES mpfr.h)

    # Restore original suffixes so other things don't break
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

    if(NOT (GMP_LIB AND GMP_INCLUDE_DIR))
        message(FATAL_ERROR "GMP static lib not found.")
    endif()

    if(NOT (MPFR_LIB AND MPFR_INCLUDE_DIR))
        message(FATAL_ERROR "MPFR static lib not found.")
    endif()
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
        message(FATAL_ERROR "FLINT not found. Ensure 'vcpkg install flint:x64-windows-static' ran.")
    endif()
else()
    # LINUX / MACOS: Build FLINT from Source (Tarball)
    # This ensures we get a static .a file linked into our binary
    include(ExternalProject)
    
    get_filename_component(GMP_LIB_DIR ${GMP_LIB} DIRECTORY)
    get_filename_component(MPFR_LIB_DIR ${MPFR_LIB} DIRECTORY)

    ExternalProject_Add(flint_build
        URL https://github.com/flintlib/flint/releases/download/v3.0.1/flint-3.0.1.tar.gz
        PREFIX ${CMAKE_BINARY_DIR}/third_party/flint
        BUILD_IN_SOURCE 1
        # Configure FLINT to use the system's (static) GMP/MPFR
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

if(NOT WIN32)
    add_dependencies(fracessa_lib flint_build)
    find_package(Threads REQUIRED)
    target_link_libraries(fracessa_lib PUBLIC Threads::Threads)
endif()

add_executable(fracessa src/main.cpp)

target_link_libraries(fracessa PRIVATE
    fracessa_lib
    argparse::argparse
)

install(TARGETS fracessa RUNTIME DESTINATION bin)
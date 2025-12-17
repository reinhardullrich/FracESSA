cmake_minimum_required(VERSION 3.10.2)

# FORCE Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Maximum optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "C++ Release flags" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "C Release flags" FORCE)

project(fracessa
    VERSION 1.0.0
    DESCRIPTION "FRACESSA - Fractional ESS Analyzer"
    LANGUAGES CXX C
)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release, but allow Debug builds
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Set Debug flags if Debug build is requested
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "C++ Debug flags" FORCE)
    set(CMAKE_C_FLAGS_DEBUG "-g -O0" CACHE STRING "C Debug flags" FORCE)
    message(STATUS "Building in Debug mode")
endif()

# LTO if supported
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND NOT MSVC)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)
    if(IPO_SUPPORTED)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
        message(STATUS "LTO/IPO enabled")
    endif()
endif()

# Compiler cache
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    # Optionally set ccache cache directory to project .cache folder
    get_filename_component(PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
    set(CCACHE_DIR "${PROJECT_ROOT}/.cache/ccache")
    if(EXISTS "${CCACHE_DIR}")
        set(ENV{CCACHE_DIR} "${CCACHE_DIR}")
        message(STATUS "ccache found, using cache directory: ${CCACHE_DIR}")
    else()
        message(STATUS "ccache found, using default cache directory")
    endif()
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Parallel builds
include(ProcessorCount)
ProcessorCount(N)
if(NOT N EQUAL 0)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N})
endif()

# Options
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_TESTING "Enable testing" OFF)
option(BUILD_TESTING "Build unit tests" OFF)
option(ENABLE_STATIC_LINKING "Enable static linking" OFF)

# Warnings
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# FetchContent
set(FETCHCONTENT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
include(FetchContent)

# Dependencies
FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.13.0
)

FetchContent_Declare(argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v2.9
)

FetchContent_MakeAvailable(
    spdlog argparse
)

# FLINT requires autotools build on Linux (CMake is Windows-only)
# Use ExternalProject to build FLINT with autotools
# FLINT requires both GMP and MPFR
include(ExternalProject)

# ============================================================================
# MPFR Detection (FLINT dependency) - Windows, Linux, macOS
# ============================================================================

# Initialize
set(MPFR_FOUND FALSE)
set(MPFR_INCLUDE_DIR "")
set(MPFR_LIB "")

# Platform-specific search
if(WIN32 OR MSVC OR MINGW OR MSYS)
    # Windows: Check standard locations
    find_path(MPFR_INCLUDE_DIR NAMES mpfr.h
        PATHS
            "C:/msys64/mingw64/include"
            "C:/gmp-install/include"
            "C:/Program Files/mpfr/include"
    )
    find_library(MPFR_LIB NAMES mpfr libmpfr
        PATHS
            "C:/msys64/mingw64/lib"
            "C:/gmp-install/lib"
            "C:/Program Files/mpfr/lib"
    )
elseif(APPLE)
    # macOS: Try pkg-config first, then Homebrew paths
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(MPFR QUIET mpfr)
        if(MPFR_FOUND AND MPFR_INCLUDE_DIRS)
            list(GET MPFR_INCLUDE_DIRS 0 MPFR_INCLUDE_DIR)
        endif()
    endif()
    
    if(NOT MPFR_INCLUDE_DIR)
        find_path(MPFR_INCLUDE_DIR NAMES mpfr.h
            PATHS
                /opt/homebrew/include
                /usr/local/include
                /opt/local/include
        )
    endif()
    
    find_library(MPFR_LIB NAMES mpfr
        PATHS
            /opt/homebrew/lib
            /usr/local/lib
            /opt/local/lib
    )
else()
    # Linux: Try pkg-config first, then standard locations
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(MPFR QUIET mpfr)
        if(MPFR_FOUND AND MPFR_INCLUDE_DIRS)
            list(GET MPFR_INCLUDE_DIRS 0 MPFR_INCLUDE_DIR)
        endif()
    endif()
    
    # Always do manual search to ensure we find both include and library
    if(NOT MPFR_INCLUDE_DIR)
        # Try common locations directly first
        if(EXISTS "/usr/include/mpfr.h")
            set(MPFR_INCLUDE_DIR "/usr/include")
        elseif(EXISTS "/usr/include/x86_64-linux-gnu/mpfr.h")
            set(MPFR_INCLUDE_DIR "/usr/include/x86_64-linux-gnu")
        else()
            find_path(MPFR_INCLUDE_DIR NAMES mpfr.h
                PATHS
                    /usr/include
                    /usr/include/x86_64-linux-gnu
                    /usr/local/include
            )
        endif()
    endif()
    
    # Always search for library manually
    if(NOT MPFR_LIB)
        # If pkg-config provided libraries, try to extract the path
        if(MPFR_FOUND AND MPFR_LIBRARIES)
            # pkg-config might return -lmpfr, try to find the actual library
            find_library(MPFR_LIB NAMES mpfr
                PATHS
                    /usr/lib
                    /usr/lib/x86_64-linux-gnu
                    /usr/lib64
                    /usr/local/lib
            )
        endif()
        
        # If still not found, try common locations directly
        if(NOT MPFR_LIB)
            if(EXISTS "/usr/lib/x86_64-linux-gnu/libmpfr.so")
                set(MPFR_LIB "/usr/lib/x86_64-linux-gnu/libmpfr.so")
            elseif(EXISTS "/usr/lib/libmpfr.so")
                set(MPFR_LIB "/usr/lib/libmpfr.so")
            else()
                find_library(MPFR_LIB NAMES mpfr
                    PATHS
                        /usr/lib
                        /usr/lib/x86_64-linux-gnu
                        /usr/lib64
                        /usr/local/lib
                )
            endif()
        endif()
    endif()
endif()

if(MPFR_INCLUDE_DIR AND MPFR_LIB)
    set(MPFR_FOUND TRUE)
    message(STATUS "Found MPFR:")
    message(STATUS "  Include: ${MPFR_INCLUDE_DIR}")
    message(STATUS "  Library: ${MPFR_LIB}")
else()
    set(MPFR_FOUND FALSE)
    message(STATUS "MPFR not found. FLINT requires MPFR.")
    message(STATUS "  Linux:   sudo apt-get install libmpfr-dev")
    message(STATUS "  macOS:   brew install mpfr")
    message(STATUS "  Windows: Install via MSYS2 or vcpkg")
endif()

# ============================================================================
# GMP Detection - Windows, Linux, macOS
# ============================================================================
# Priority: System GMP first, then local GMP as fallback

set(GMP_FOUND FALSE)
set(GMP_INCLUDE_DIR "")
set(GMP_LIBRARIES "")

# Platform-specific search (system GMP first)
if(WIN32 OR MSVC OR MINGW OR MSYS)
    # ===== WINDOWS =====
    message(STATUS "Detecting GMP on Windows...")
    
    # Check for manually built GMP first (from MSYS2/MinGW build)
    # Environment variables set by workflow
    if(DEFINED ENV{GMP_ROOT})
        set(GMP_ROOT $ENV{GMP_ROOT})
        list(APPEND GMP_SEARCH_PATHS "${GMP_ROOT}/include")
        list(APPEND GMP_LIB_SEARCH_PATHS "${GMP_ROOT}/lib")
        message(STATUS "Using GMP_ROOT from environment: ${GMP_ROOT}")
    endif()
    
    if(DEFINED ENV{GMP_INCLUDE_DIR})
        list(INSERT GMP_SEARCH_PATHS 0 "$ENV{GMP_INCLUDE_DIR}")
        message(STATUS "Using GMP_INCLUDE_DIR from environment: $ENV{GMP_INCLUDE_DIR}")
    endif()
    
    if(DEFINED ENV{GMP_LIBRARY_DIR})
        list(INSERT GMP_LIB_SEARCH_PATHS 0 "$ENV{GMP_LIBRARY_DIR}")
        message(STATUS "Using GMP_LIBRARY_DIR from environment: $ENV{GMP_LIBRARY_DIR}")
    endif()
    
    # Add standard GMP locations (MSYS2 package installs to /mingw64)
    # Prioritize MSYS2/MinGW paths where pacman installs GMP
    list(INSERT GMP_SEARCH_PATHS 0
        "C:/msys64/mingw64/include"
    )
    list(INSERT GMP_LIB_SEARCH_PATHS 0
        "C:/msys64/mingw64/lib"
    )
    # Add other fallback paths
    list(APPEND GMP_SEARCH_PATHS
        "C:/gmp-install/include"
        "C:/Program Files/gmp/include"
        "C:/gmp/include"
    )
    list(APPEND GMP_LIB_SEARCH_PATHS
        "C:/gmp-install/lib"
        "C:/Program Files/gmp/lib"
        "C:/gmp/lib"
    )
    
    # Check vcpkg as fallback (if available)
    if(DEFINED ENV{VCPKG_ROOT})
        set(VCPKG_ROOT $ENV{VCPKG_ROOT})
        message(STATUS "VCPKG_ROOT: ${VCPKG_ROOT} (checking as fallback)")
        
        # Use VCPKG_TARGET_TRIPLET from environment if set
        if(DEFINED ENV{VCPKG_TARGET_TRIPLET})
            set(VCPKG_TRIPLET $ENV{VCPKG_TARGET_TRIPLET})
        else()
            if(CMAKE_SIZEOF_VOID_P EQUAL 8)
                set(VCPKG_TRIPLET "x64-windows")
            else()
                set(VCPKG_TRIPLET "x86-windows")
            endif()
            if(ENABLE_STATIC_LINKING)
                set(VCPKG_TRIPLET "${VCPKG_TRIPLET}-static")
            endif()
        endif()
        
        # Add vcpkg paths as fallback
        list(APPEND GMP_SEARCH_PATHS
            "${VCPKG_ROOT}/installed/${VCPKG_TRIPLET}/include"
            "${VCPKG_ROOT}/vcpkg_installed/${VCPKG_TRIPLET}/include"
        )
        list(APPEND GMP_LIB_SEARCH_PATHS
            "${VCPKG_ROOT}/installed/${VCPKG_TRIPLET}/lib"
            "${VCPKG_ROOT}/vcpkg_installed/${VCPKG_TRIPLET}/lib"
        )
    endif()
    
    # Check if GMP variables are already set (e.g., from command line -D flags)
    if(GMP_INCLUDE_DIR AND GMP_LIB)
        message(STATUS "Using GMP variables set from command line:")
        message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
        message(STATUS "  Library: ${GMP_LIB}")
        set(GMP_FOUND TRUE)
        if(NOT GMP_LIBRARIES)
            set(GMP_LIBRARIES ${GMP_LIB})
        endif()
    else()
        message(STATUS "Searching for GMP in:")
        message(STATUS "  Include paths: ${GMP_SEARCH_PATHS}")
        message(STATUS "  Library paths: ${GMP_LIB_SEARCH_PATHS}")
        
        # Search for headers (gmp.h - only C library needed for FLINT)
        find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS ${GMP_SEARCH_PATHS} NO_DEFAULT_PATH)
        if(NOT GMP_INCLUDE_DIR)
            find_path(GMP_INCLUDE_DIR NAMES gmp.h) # Fallback
        endif()
        
        # Search for libraries (Windows/MinGW uses .dll.a files for dynamic linking)
        # Prioritize .dll.a (import libraries) for dynamic linking, fallback to .a (static)
        find_library(GMP_LIB 
            NAMES libgmp.dll.a gmp.dll.a libgmp gmp libgmp.a gmp.a gmp.lib libgmp.lib
            PATHS ${GMP_LIB_SEARCH_PATHS} 
            NO_DEFAULT_PATH
        )
        if(NOT GMP_LIB)
            find_library(GMP_LIB NAMES libgmp.dll.a gmp.dll.a libgmp gmp libgmp.a gmp.a gmp.lib libgmp.lib)
        endif()
        
        if(GMP_INCLUDE_DIR AND GMP_LIB)
            set(GMP_FOUND TRUE)
            set(GMP_LIBRARIES ${GMP_LIB})
            message(STATUS "Found GMP on Windows:")
            message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
            message(STATUS "  Library: ${GMP_LIB}")
        else()
            message(STATUS "GMP detection failed:")
            message(STATUS "  Include dir: ${GMP_INCLUDE_DIR}")
            message(STATUS "  Library: ${GMP_LIB}")
        endif()
    endif()
    
elseif(APPLE)
    # ===== macOS =====
    message(STATUS "Detecting GMP on macOS...")
    
    # Try pkg-config first (Homebrew often provides .pc files)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(GMP QUIET gmp)
        if(GMP_FOUND AND GMP_INCLUDE_DIRS)
            list(GET GMP_INCLUDE_DIRS 0 GMP_INCLUDE_DIR)
            message(STATUS "Found GMP via pkg-config: ${GMP_INCLUDE_DIR}")
        endif()
    endif()
    
    # Build search paths (environment variables take precedence)
    set(GMP_SEARCH_PATHS "")
    set(GMP_LIB_SEARCH_PATHS "")
    
    # Check environment variables first (set by workflow)
    if(DEFINED ENV{GMP_INCLUDE_DIR})
        list(APPEND GMP_SEARCH_PATHS "$ENV{GMP_INCLUDE_DIR}")
        message(STATUS "Using GMP_INCLUDE_DIR from environment: $ENV{GMP_INCLUDE_DIR}")
    elseif(DEFINED ENV{GMP_ROOT})
        list(APPEND GMP_SEARCH_PATHS "$ENV{GMP_ROOT}/include")
        message(STATUS "Using GMP_ROOT from environment: $ENV{GMP_ROOT}")
    endif()
    
    if(DEFINED ENV{GMP_LIBRARY_DIR})
        list(APPEND GMP_LIB_SEARCH_PATHS "$ENV{GMP_LIBRARY_DIR}")
        message(STATUS "Using GMP_LIBRARY_DIR from environment: $ENV{GMP_LIBRARY_DIR}")
    elseif(DEFINED ENV{GMP_ROOT})
        list(APPEND GMP_LIB_SEARCH_PATHS "$ENV{GMP_ROOT}/lib")
    endif()
    
    # Add standard Homebrew paths (check both Apple Silicon and Intel)
    list(APPEND GMP_SEARCH_PATHS
        /opt/homebrew/include  # Apple Silicon
        /usr/local/include     # Intel
        /opt/local/include     # MacPorts
    )
    list(APPEND GMP_LIB_SEARCH_PATHS
        /opt/homebrew/lib
        /usr/local/lib
        /opt/local/lib
    )
    
    # Only search if not found via pkg-config
    if(NOT GMP_INCLUDE_DIR)
        find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS ${GMP_SEARCH_PATHS} NO_DEFAULT_PATH)
        if(NOT GMP_INCLUDE_DIR)
            find_path(GMP_INCLUDE_DIR NAMES gmp.h)
        endif()
    endif()
    
    # Always search for libraries manually (pkg-config may not provide full paths)
    find_library(GMP_LIB NAMES gmp PATHS ${GMP_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(NOT GMP_LIB)
        find_library(GMP_LIB NAMES gmp)
    endif()
    
    if(GMP_INCLUDE_DIR AND GMP_LIB)
        set(GMP_FOUND TRUE)
        set(GMP_LIBRARIES ${GMP_LIB})
        message(STATUS "Found GMP on macOS:")
        message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
        message(STATUS "  GMP lib: ${GMP_LIB}")
    else()
        message(STATUS "GMP detection failed:")
        message(STATUS "  Include dir: ${GMP_INCLUDE_DIR}")
        message(STATUS "  GMP lib: ${GMP_LIB}")
    endif()
    
else()
    # ===== LINUX =====
    message(STATUS "Detecting GMP on Linux...")
    
    # Try pkg-config first for include directory
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(GMP QUIET gmp)
        if(GMP_FOUND AND GMP_INCLUDE_DIRS)
            list(GET GMP_INCLUDE_DIRS 0 GMP_INCLUDE_DIR)
        endif()
    endif()
    
    # Manual search for include directory if not found via pkg-config
    if(NOT GMP_INCLUDE_DIR)
        # Try common locations directly first
        if(EXISTS "/usr/include/gmp.h")
            set(GMP_INCLUDE_DIR "/usr/include")
        elseif(EXISTS "/usr/include/x86_64-linux-gnu/gmp.h")
            set(GMP_INCLUDE_DIR "/usr/include/x86_64-linux-gnu")
        else()
            set(GMP_SEARCH_PATHS
                /usr/include
                /usr/include/x86_64-linux-gnu
                /usr/local/include
            )
            # Find gmp.h (C library only, needed for FLINT)
            find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS ${GMP_SEARCH_PATHS} NO_DEFAULT_PATH)
            if(NOT GMP_INCLUDE_DIR)
                find_path(GMP_INCLUDE_DIR NAMES gmp.h)
            endif()
        endif()
    endif()
    
    # Always search for libraries manually
    set(GMP_LIB_SEARCH_PATHS
        /usr/lib
        /usr/lib/x86_64-linux-gnu
        /usr/lib64
        /usr/local/lib
    )
    
    find_library(GMP_LIB NAMES gmp PATHS ${GMP_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(NOT GMP_LIB)
        find_library(GMP_LIB NAMES gmp)
    endif()
    
    if(GMP_INCLUDE_DIR AND GMP_LIB)
        set(GMP_FOUND TRUE)
        set(GMP_LIBRARIES ${GMP_LIB})
        if(PkgConfig_FOUND)
            message(STATUS "Found GMP via pkg-config and manual library search")
        else()
            message(STATUS "Found GMP manually on Linux")
        endif()
    else()
        set(GMP_FOUND FALSE)
        message(STATUS "GMP detection failed - missing components")
    endif()
endif()  # End of platform-specific search

# Fallback to local GMP if system GMP not found
if(NOT GMP_FOUND)
    get_filename_component(PROJECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
    set(PROJECT_GMP_DIR "${PROJECT_ROOT}/gmp")
    if(EXISTS "${PROJECT_GMP_DIR}/include/gmp.h")
        set(GMP_INCLUDE_DIR "${PROJECT_GMP_DIR}/include")
        find_library(GMP_LIB NAMES gmp PATHS ${PROJECT_GMP_DIR}/lib NO_DEFAULT_PATH)
        if(GMP_LIB)
            set(GMP_FOUND TRUE)
            set(GMP_LIBRARIES ${GMP_LIB})
            message(STATUS "Found local GMP 6.3 in project directory (fallback):")
            message(STATUS "  Path: ${PROJECT_GMP_DIR}")
            message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
            message(STATUS "  Libraries: ${GMP_LIBRARIES}")
        endif()
    endif()
endif()

# Verify GMP found
if(NOT GMP_FOUND OR NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES)
    message(FATAL_ERROR "
    GMP library not found!
    
    Installation instructions:
    - Linux:   sudo apt-get install libgmp-dev
    - macOS:   brew install gmp
    - Windows: vcpkg install gmp
               (or download from https://gmplib.org/)
    ")
endif()

message(STATUS "GMP Configuration:")
message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
message(STATUS "  Libraries: ${GMP_LIBRARIES}")

# Verify header (check for gmp.h)
if(NOT EXISTS "${GMP_INCLUDE_DIR}/gmp.h")
    message(FATAL_ERROR "GMP header not found in: ${GMP_INCLUDE_DIR}")
endif()

# Static linking
if(ENABLE_STATIC_LINKING)
    message(STATUS "Static linking requested")
    
    if(WIN32)
        # Windows: vcpkg handles this automatically with triplet
        message(STATUS "Windows: Static linking via vcpkg triplet")
    else()
        # Linux/macOS: find static libraries
        get_filename_component(GMP_LIB_DIR "${GMP_LIB}" DIRECTORY)
        find_library(GMP_STATIC NAMES libgmp.a PATHS ${GMP_LIB_DIR} NO_DEFAULT_PATH)
        if(GMP_STATIC)
            set(GMP_LIBRARIES ${GMP_STATIC})
            message(STATUS "Using static GMP")
            
            if(UNIX AND NOT APPLE)
                # Static libgcc/libstdc++ but dynamic glibc
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
                message(STATUS "Static linking: libgcc, libstdc++ (glibc dynamic for compatibility)")
            endif()
        else()
            message(WARNING "Static GMP not found, using dynamic")
        endif()
    endif()
endif()

# ============================================================================
# FLINT Build Configuration
# ============================================================================
# FLINT requires both GMP and MPFR

set(FLINT_BUILD_ENABLED FALSE)

if(GMP_FOUND AND MPFR_FOUND)
    set(FLINT_BUILD_ENABLED TRUE)
    
    # Determine GMP root directory for FLINT configure
    get_filename_component(GMP_LIB_DIR "${GMP_LIB}" DIRECTORY)
    get_filename_component(GMP_ROOT_DIR "${GMP_LIB_DIR}/.." ABSOLUTE)
    
    # Determine MPFR root directory for FLINT configure
    get_filename_component(MPFR_LIB_DIR "${MPFR_LIB}" DIRECTORY)
    get_filename_component(MPFR_ROOT_DIR "${MPFR_LIB_DIR}/.." ABSOLUTE)
    
    # Build FLINT using ExternalProject with autotools
    # Using v3.0.0 (compatible with system GMP 6.2.0 and MPFR 4.0.2) with both GMP and MPFR
    ExternalProject_Add(flint_build
        GIT_REPOSITORY https://github.com/flintlib/flint.git
        GIT_TAG v3.0.0
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/flint-src
        CONFIGURE_COMMAND bash -c "cd <SOURCE_DIR> && ./bootstrap.sh && CPPFLAGS='-I${GMP_INCLUDE_DIR} -I${MPFR_INCLUDE_DIR}' LDFLAGS='-L${GMP_LIB_DIR} -L${MPFR_LIB_DIR}' ./configure --prefix=<INSTALL_DIR> --with-gmp=${GMP_ROOT_DIR} --with-mpfr=${MPFR_ROOT_DIR} --enable-shared=no --enable-static=yes"
        BUILD_COMMAND bash -c "cd <SOURCE_DIR> && make -j$(nproc)"
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND bash -c "cd <SOURCE_DIR> && make install"
        BUILD_ALWAYS FALSE
        UPDATE_COMMAND ""
    )
    
    # Set FLINT variables for linking
    ExternalProject_Get_Property(flint_build INSTALL_DIR)
    set(FLINT_LIB_DIR ${INSTALL_DIR}/lib)
    set(FLINT_INCLUDE_DIR ${INSTALL_DIR}/include)
    set(FLINT_LIBRARIES ${FLINT_LIB_DIR}/libflint.a)
    
    message(STATUS "FLINT will be built with GMP and MPFR")
    message(STATUS "FLINT will be installed to: ${INSTALL_DIR}")
else()
    if(NOT GMP_FOUND)
        message(STATUS "FLINT build disabled: GMP not found")
    endif()
    if(NOT MPFR_FOUND)
        message(STATUS "FLINT build disabled: MPFR not found")
    endif()
endif()

# ============================================================================
# Build Library
# ============================================================================

add_library(fracessa_lib
    src/fracessa.cpp
    src/findeq.cpp
    src/checkstab.cpp
)

set_target_properties(fracessa_lib PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(fracessa_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/fracessa>
        $<BUILD_INTERFACE:${argparse_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${spdlog_SOURCE_DIR}/include>
        ${GMP_INCLUDE_DIR}
        ${MPFR_INCLUDE_DIR}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Add FLINT include directory if enabled
if(FLINT_BUILD_ENABLED)
    target_include_directories(fracessa_lib PRIVATE
        ${FLINT_INCLUDE_DIR}
    )
    # Add dependency so FLINT builds before fracessa_lib
    add_dependencies(fracessa_lib flint_build)
endif()

target_link_libraries(fracessa_lib PRIVATE
    ${GMP_LIBRARIES}
    ${MPFR_LIB}
)

# Link FLINT library if enabled
if(FLINT_BUILD_ENABLED)
    target_link_libraries(fracessa_lib PRIVATE
        ${FLINT_LIBRARIES}
    )
endif()

# ============================================================================
# Build Executable
# ============================================================================

add_executable(fracessa src/main.cpp)

target_include_directories(fracessa PRIVATE
    ${argparse_SOURCE_DIR}/include
    ${GMP_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
)

# Add FLINT include directory if enabled
if(FLINT_BUILD_ENABLED)
    target_include_directories(fracessa PRIVATE
        ${FLINT_INCLUDE_DIR}
    )
endif()

target_link_libraries(fracessa PRIVATE
    fracessa_lib
    ${GMP_LIBRARIES}
    ${MPFR_LIB}
)

# Add FLINT library and pthread if enabled
if(FLINT_BUILD_ENABLED)
    target_link_libraries(fracessa PRIVATE
        ${FLINT_LIBRARIES}
        pthread
    )
endif()

# Windows-specific
if(WIN32 OR MSVC OR MSYS OR MINGW)
    if(ENABLE_STATIC_LINKING)
        set_target_properties(fracessa PROPERTIES
            LINK_FLAGS "-static"
        )
    endif()
endif()

# Install
install(TARGETS fracessa fracessa_lib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/fracessa/
    DESTINATION include/fracessa
    FILES_MATCHING PATTERN "*.hpp"
)

# Testing
if(BUILD_TESTING)
    enable_testing()
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
endif()
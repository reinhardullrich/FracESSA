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

# Enforce Release
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
    message(FATAL_ERROR "Only Release builds supported. Current: ${CMAKE_BUILD_TYPE}")
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

# ============================================================================
# GMP Detection - Windows, Linux, macOS
# ============================================================================

# Initialize
set(GMP_FOUND FALSE)
set(GMP_INCLUDE_DIR "")
set(GMP_LIBRARIES "")

# Platform-specific search
if(WIN32 OR MSVC OR MINGW OR MSYS)
    # ===== WINDOWS =====
    message(STATUS "Detecting GMP on Windows...")
    
    # Check vcpkg (most common on Windows)
    if(DEFINED ENV{VCPKG_ROOT})
        set(VCPKG_ROOT $ENV{VCPKG_ROOT})
        message(STATUS "VCPKG_ROOT: ${VCPKG_ROOT}")
        
        # Use VCPKG_TARGET_TRIPLET from environment if set (set by vcpkg action)
        if(DEFINED ENV{VCPKG_TARGET_TRIPLET})
            set(VCPKG_TRIPLET $ENV{VCPKG_TARGET_TRIPLET})
            message(STATUS "Using VCPKG_TARGET_TRIPLET from environment: ${VCPKG_TRIPLET}")
        else()
            # Auto-detect triplet as fallback
            if(CMAKE_SIZEOF_VOID_P EQUAL 8)
                set(VCPKG_TRIPLET "x64-windows")
            else()
                set(VCPKG_TRIPLET "x86-windows")
            endif()
            
            if(ENABLE_STATIC_LINKING)
                set(VCPKG_TRIPLET "${VCPKG_TRIPLET}-static")
            endif()
            message(STATUS "Auto-detected vcpkg triplet: ${VCPKG_TRIPLET}")
        endif()
        
        # Search in both installed and vcpkg_installed directories
        # vcpkg_installed is used by newer vcpkg versions
        set(GMP_SEARCH_PATHS
            "${VCPKG_ROOT}/installed/${VCPKG_TRIPLET}/include"
            "${VCPKG_ROOT}/vcpkg_installed/${VCPKG_TRIPLET}/include"
            "${VCPKG_ROOT}/packages/gmp_${VCPKG_TRIPLET}/include"
        )
        
        set(GMP_LIB_SEARCH_PATHS
            "${VCPKG_ROOT}/installed/${VCPKG_TRIPLET}/lib"
            "${VCPKG_ROOT}/installed/${VCPKG_TRIPLET}/debug/lib"
            "${VCPKG_ROOT}/vcpkg_installed/${VCPKG_TRIPLET}/lib"
            "${VCPKG_ROOT}/vcpkg_installed/${VCPKG_TRIPLET}/debug/lib"
            "${VCPKG_ROOT}/packages/gmp_${VCPKG_TRIPLET}/lib"
        )
        
        message(STATUS "Searching for GMP in:")
        message(STATUS "  Include paths: ${GMP_SEARCH_PATHS}")
        message(STATUS "  Library paths: ${GMP_LIB_SEARCH_PATHS}")
    else()
        # Fallback Windows paths
        set(GMP_SEARCH_PATHS
            "C:/Program Files/gmp/include"
            "C:/gmp/include"
            "C:/msys64/mingw64/include"
        )
        set(GMP_LIB_SEARCH_PATHS
            "C:/Program Files/gmp/lib"
            "C:/gmp/lib"
            "C:/msys64/mingw64/lib"
        )
    endif()
    
    # Search for headers (try both gmp.h and gmpxx.h)
    find_path(GMP_INCLUDE_DIR NAMES gmpxx.h gmp.h PATHS ${GMP_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(NOT GMP_INCLUDE_DIR)
        find_path(GMP_INCLUDE_DIR NAMES gmpxx.h gmp.h) # Fallback
    endif()
    
    # Search for libraries (Windows uses .lib files)
    find_library(GMP_LIB 
        NAMES gmp libgmp gmp.lib libgmp.lib
        PATHS ${GMP_LIB_SEARCH_PATHS} 
        NO_DEFAULT_PATH
    )
    if(NOT GMP_LIB)
        find_library(GMP_LIB NAMES gmp libgmp gmp.lib libgmp.lib)
    endif()
    
    find_library(GMPXX_LIB 
        NAMES gmpxx libgmpxx gmpxx.lib libgmpxx.lib
        PATHS ${GMP_LIB_SEARCH_PATHS} 
        NO_DEFAULT_PATH
    )
    if(NOT GMPXX_LIB)
        find_library(GMPXX_LIB NAMES gmpxx libgmpxx gmpxx.lib libgmpxx.lib)
    endif()
    
    if(GMP_INCLUDE_DIR AND GMP_LIB AND GMPXX_LIB)
        set(GMP_FOUND TRUE)
        set(GMP_LIBRARIES ${GMP_LIB} ${GMPXX_LIB})
        message(STATUS "Found GMP on Windows:")
        message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
        message(STATUS "  GMP lib: ${GMP_LIB}")
        message(STATUS "  GMPXX lib: ${GMPXX_LIB}")
    else()
        message(STATUS "GMP detection failed:")
        message(STATUS "  Include dir: ${GMP_INCLUDE_DIR}")
        message(STATUS "  GMP lib: ${GMP_LIB}")
        message(STATUS "  GMPXX lib: ${GMPXX_LIB}")
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
        find_path(GMP_INCLUDE_DIR NAMES gmpxx.h gmp.h PATHS ${GMP_SEARCH_PATHS} NO_DEFAULT_PATH)
        if(NOT GMP_INCLUDE_DIR)
            find_path(GMP_INCLUDE_DIR NAMES gmpxx.h gmp.h)
        endif()
    endif()
    
    # Always search for libraries manually (pkg-config may not provide full paths)
    find_library(GMP_LIB NAMES gmp PATHS ${GMP_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(NOT GMP_LIB)
        find_library(GMP_LIB NAMES gmp)
    endif()
    
    find_library(GMPXX_LIB NAMES gmpxx PATHS ${GMP_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(NOT GMPXX_LIB)
        find_library(GMPXX_LIB NAMES gmpxx)
    endif()
    
    if(GMP_INCLUDE_DIR AND GMP_LIB AND GMPXX_LIB)
        set(GMP_FOUND TRUE)
        set(GMP_LIBRARIES ${GMP_LIB} ${GMPXX_LIB})
        message(STATUS "Found GMP on macOS:")
        message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
        message(STATUS "  GMP lib: ${GMP_LIB}")
        message(STATUS "  GMPXX lib: ${GMPXX_LIB}")
    else()
        message(STATUS "GMP detection failed:")
        message(STATUS "  Include dir: ${GMP_INCLUDE_DIR}")
        message(STATUS "  GMP lib: ${GMP_LIB}")
        message(STATUS "  GMPXX lib: ${GMPXX_LIB}")
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
        if(EXISTS "/usr/include/gmpxx.h")
            set(GMP_INCLUDE_DIR "/usr/include")
        elseif(EXISTS "/usr/include/x86_64-linux-gnu/gmpxx.h")
            set(GMP_INCLUDE_DIR "/usr/include/x86_64-linux-gnu")
        else()
            set(GMP_SEARCH_PATHS
                /usr/include
                /usr/include/x86_64-linux-gnu
                /usr/local/include
            )
            # Try to find gmpxx.h first, then fall back to gmp.h
            find_path(GMP_INCLUDE_DIR NAMES gmpxx.h PATHS ${GMP_SEARCH_PATHS} NO_DEFAULT_PATH)
            if(NOT GMP_INCLUDE_DIR)
                find_path(GMP_INCLUDE_DIR NAMES gmpxx.h)
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
    find_library(GMPXX_LIB NAMES gmpxx PATHS ${GMP_LIB_SEARCH_PATHS} NO_DEFAULT_PATH)
    if(NOT GMPXX_LIB)
        find_library(GMPXX_LIB NAMES gmpxx)
    endif()
    
    if(GMP_INCLUDE_DIR AND GMP_LIB AND GMPXX_LIB)
        set(GMP_FOUND TRUE)
        set(GMP_LIBRARIES ${GMP_LIB} ${GMPXX_LIB})
        if(PkgConfig_FOUND)
            message(STATUS "Found GMP via pkg-config and manual library search")
        else()
            message(STATUS "Found GMP manually on Linux")
        endif()
    else()
        set(GMP_FOUND FALSE)
        message(STATUS "GMP detection failed - missing components")
    endif()
endif()

# Verify GMP found
if(NOT GMP_FOUND OR NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES)
    message(FATAL_ERROR "
    GMP library not found!
    
    Installation instructions:
    - Linux:   sudo apt-get install libgmp-dev libgmpxx4ldbl
    - macOS:   brew install gmp
    - Windows: vcpkg install gmp gmpxx
               (or download from https://gmplib.org/)
    ")
endif()

message(STATUS "GMP Configuration:")
message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
message(STATUS "  Libraries: ${GMP_LIBRARIES}")

# Verify header (check for gmpxx.h which includes gmp.h)
if(NOT EXISTS "${GMP_INCLUDE_DIR}/gmpxx.h" AND NOT EXISTS "${GMP_INCLUDE_DIR}/gmp.h")
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
        find_library(GMPXX_STATIC NAMES libgmpxx.a PATHS ${GMP_LIB_DIR} NO_DEFAULT_PATH)
        
        if(GMP_STATIC AND GMPXX_STATIC)
            set(GMP_LIBRARIES ${GMP_STATIC} ${GMPXX_STATIC})
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
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)


target_link_libraries(fracessa_lib PRIVATE
    ${GMP_LIBRARIES}
)

# ============================================================================
# Build Executable
# ============================================================================

add_executable(fracessa src/main.cpp)

target_include_directories(fracessa PRIVATE
    ${argparse_SOURCE_DIR}/include
    ${GMP_INCLUDE_DIR}
)

target_link_libraries(fracessa PRIVATE
    fracessa_lib
    ${GMP_LIBRARIES}
)

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
cmake_minimum_required(VERSION 3.18)

project(exact_integer_solver_comparison LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_path(FLINT_INCLUDE_DIR NAMES flint/fmpq_mat.h REQUIRED)
find_library(FLINT_LIBRARY NAMES flint REQUIRED)
find_library(MPFR_LIBRARY NAMES mpfr REQUIRED)
find_library(GMP_LIBRARY NAMES gmp REQUIRED)

add_executable(exact_integer_solver_benchmark
    benchmark.cpp
    ../../cpp/src/find_candidate_exact.cpp
)

target_include_directories(exact_integer_solver_benchmark PRIVATE
    ../../cpp/include
    ${FLINT_INCLUDE_DIR}
)

target_link_libraries(exact_integer_solver_benchmark PRIVATE
    ${FLINT_LIBRARY}
    ${MPFR_LIBRARY}
    ${GMP_LIBRARY}
)

include(CTest)
if(BUILD_TESTING)
    # The full-support reduced system is [[0,1],[1,0]]: nonsingular and indefinite. FFLDLT must return its candidate for superset pruning, not reject it.
    add_test(NAME FFLDLTIndefiniteCandidate COMMAND exact_integer_solver_benchmark "3#0,-1/4,-1/4,-1/2,1/2,-1/2" 1)
endif()

if(NOT MSVC)
    target_compile_options(exact_integer_solver_benchmark PRIVATE
        $<$<CONFIG:Release>:-funroll-loops>
        $<$<CONFIG:Release>:-march=native>
    )

    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED)
    if(IPO_SUPPORTED)
        set_property(TARGET exact_integer_solver_benchmark PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
endif()

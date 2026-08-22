cmake_minimum_required(VERSION 3.18)

project(fracessa_a1 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(FRACESSA_BUILD_CLI ON CACHE BOOL "" FORCE)
set(FRACESSA_BUILD_PYTHON ON CACHE BOOL "" FORCE)
find_package(Python3 3.11 COMPONENTS Interpreter Development.Module REQUIRED)

# Build an unchanged production baseline inside this model's build tree.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../cpp" "${CMAKE_CURRENT_BINARY_DIR}/baseline" EXCLUDE_FROM_ALL)

add_library(fracessa_a1_lib
    src/fracessa.cpp
    ../../cpp/src/analyzer_logging.cpp
    src/fast_candidate_filter.cpp
    src/exact_candidate_solver.cpp
    ../../cpp/src/check_stability.cpp
)
target_include_directories(fracessa_a1_lib BEFORE PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/include"
)
target_link_libraries(fracessa_a1_lib PUBLIC coposit::safe spdlog::spdlog)

get_target_property(BASELINE_IPO fracessa_lib INTERPROCEDURAL_OPTIMIZATION_RELEASE)
if(BASELINE_IPO)
    set_property(TARGET fracessa_a1_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

if(NOT MSVC)
    target_compile_options(fracessa_a1_lib PRIVATE $<$<CONFIG:Release>:-funroll-loops>)
    if(FRACESSA_NATIVE_ARCH)
        target_compile_options(fracessa_a1_lib PRIVATE $<$<CONFIG:Release>:-march=native>)
    endif()
endif()

add_executable(fracessa_a1 ../../cpp/src/main.cpp)
target_compile_definitions(fracessa_a1 PRIVATE FRACESSA_VERSION="a1")
target_link_libraries(fracessa_a1 PRIVATE fracessa_a1_lib argparse::argparse)

pybind11_add_module(fracessa_a1_core MODULE ../../cpp/src/pybind_module.cpp)
target_link_libraries(fracessa_a1_core PRIVATE fracessa_a1_lib)
set_target_properties(fracessa_a1_core PROPERTIES
    OUTPUT_NAME fracessa_core
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/a1"
)

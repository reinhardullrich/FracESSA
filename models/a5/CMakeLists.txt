cmake_minimum_required(VERSION 3.18)

project(fracessa_a5 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(FRACESSA_BUILD_CLI ON CACHE BOOL "" FORCE)
set(FRACESSA_BUILD_PYTHON ON CACHE BOOL "" FORCE)
find_package(Python3 3.11 COMPONENTS Interpreter Development.Module REQUIRED)
include(FetchContent)

# Reuse the production dependencies and A4's unchanged numerical/adaptor sources.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../cpp" "${CMAKE_CURRENT_BINARY_DIR}/baseline" EXCLUDE_FROM_ALL)

FetchContent_Declare(cadical
    GIT_REPOSITORY https://github.com/arminbiere/cadical.git
    GIT_TAG rel-2.2.1
    GIT_SHALLOW TRUE
)
FetchContent_GetProperties(cadical)
if(NOT cadical_POPULATED)
    cmake_policy(PUSH)
    if(POLICY CMP0169)
        cmake_policy(SET CMP0169 OLD)
    endif()
    FetchContent_Populate(cadical)
    cmake_policy(POP)
endif()
file(GLOB CADICAL_LIBRARY_SOURCES CONFIGURE_DEPENDS "${cadical_SOURCE_DIR}/src/*.cpp")
list(FILTER CADICAL_LIBRARY_SOURCES EXCLUDE REGEX "/(cadical|mobical)\\.cpp$")
add_library(fracessa_a5_cadical STATIC ${CADICAL_LIBRARY_SOURCES} ${cadical_SOURCE_DIR}/src/kitten.c)
set_target_properties(fracessa_a5_cadical PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(fracessa_a5_cadical PRIVATE NBUILD)
target_include_directories(fracessa_a5_cadical PUBLIC ${cadical_SOURCE_DIR}/src)

add_library(fracessa_a5_lib
    src/fracessa.cpp
    src/exact_candidate_solver.cpp
    ../a4/src/analyzer_logging.cpp
    ../a4/src/fast_candidate_filter.cpp
    ../../cpp/src/check_stability.cpp
)
target_include_directories(fracessa_a5_lib BEFORE PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../a4/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/include"
)
target_link_libraries(fracessa_a5_lib PUBLIC coposit::safe spdlog::spdlog fracessa_a5_cadical)

get_target_property(BASELINE_IPO fracessa_lib INTERPROCEDURAL_OPTIMIZATION_RELEASE)
if(BASELINE_IPO)
    set_property(TARGET fracessa_a5_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

if(NOT MSVC)
    target_compile_options(fracessa_a5_lib PRIVATE $<$<CONFIG:Release>:-funroll-loops>)
    if(FRACESSA_NATIVE_ARCH)
        target_compile_options(fracessa_a5_lib PRIVATE $<$<CONFIG:Release>:-march=native>)
    endif()
endif()

add_executable(fracessa_a5 ../a4/src/main.cpp)
target_compile_definitions(fracessa_a5 PRIVATE FRACESSA_VERSION="a5")
target_link_libraries(fracessa_a5 PRIVATE fracessa_a5_lib argparse::argparse)

pybind11_add_module(fracessa_a5_core MODULE ../a4/src/pybind_module.cpp)
target_link_libraries(fracessa_a5_core PRIVATE fracessa_a5_lib)
set_target_properties(fracessa_a5_core PROPERTIES
    OUTPUT_NAME fracessa_core
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/a5"
)

add_executable(fracessa_a5_sat_check tests/test_support_generator_sat.cpp)
target_include_directories(fracessa_a5_sat_check PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../cpp/include"
)
target_link_libraries(fracessa_a5_sat_check PRIVATE fracessa_a5_cadical)

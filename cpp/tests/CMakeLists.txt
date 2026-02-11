# Test executables
set(TEST_LIST
    test_matrix
    test_matrix_parser
    test_supports
    test_linear_solver
    test_lu
    test_copositivity
    test_bitset64
    test_integration
    test_fraction
)

foreach(test_name ${TEST_LIST})
    add_executable(${test_name} ${test_name}.cpp)
    target_link_libraries(${test_name} PRIVATE 
        fracessa_lib 
        GTest::gtest_main 
        ${FLINT_LIBRARIES}
        ${MPFR_LIB}
        ${GMP_LIB}
    )
    target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
endforeach()

# Copy test data to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Register tests with CTest
add_test(NAME MatrixTests COMMAND test_matrix)
add_test(NAME MatrixParserTests COMMAND test_matrix_parser)
add_test(NAME SupportsTests COMMAND test_supports)
add_test(NAME LinearSolverTests COMMAND test_linear_solver)
add_test(NAME LUTests COMMAND test_lu)
add_test(NAME CopositivityTests COMMAND test_copositivity)
add_test(NAME Bitset64Tests COMMAND test_bitset64)
add_test(NAME IntegrationTests COMMAND test_integration)
add_test(NAME FractionTests COMMAND test_fraction)

# ------------------------------------------------------------------------------
# Verification correctness tests (one CTest per matrix for parallel execution)
# ------------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(CLI_PARSER_BLACKBOX_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/cli_parser_blackbox.py")
add_test(
    NAME CliParserBlackBoxTests
    COMMAND ${Python3_EXECUTABLE} ${CLI_PARSER_BLACKBOX_SCRIPT}
        --fracessa-exe $<TARGET_FILE:fracessa>
)
set_tests_properties(CliParserBlackBoxTests PROPERTIES LABELS "cli;parser")

set(VERIFICATION_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../python/verification")
set(MATRIX_SELECTOR_SCRIPT "${VERIFICATION_DIR}/matrix_selection.py")
set(VERIFY_MATRIX_SCRIPT "${VERIFICATION_DIR}/ctest_verify_matrix.py")

execute_process(
    COMMAND ${Python3_EXECUTABLE} ${MATRIX_SELECTOR_SCRIPT} --list-tests --verification-dir ${VERIFICATION_DIR}
    RESULT_VARIABLE MATRIX_LIST_RC
    OUTPUT_VARIABLE MATRIX_LIST_OUT
    ERROR_VARIABLE MATRIX_LIST_ERR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT MATRIX_LIST_RC EQUAL 0)
    message(FATAL_ERROR "Failed to get verification matrix list: ${MATRIX_LIST_ERR}")
endif()

string(REPLACE "\n" ";" MATRIX_LIST_LINES "${MATRIX_LIST_OUT}")
foreach(matrix_line IN LISTS MATRIX_LIST_LINES)
    if(matrix_line STREQUAL "")
        continue()
    endif()

    string(REPLACE "," ";" matrix_parts "${matrix_line}")
    list(LENGTH matrix_parts matrix_parts_len)
    if(NOT matrix_parts_len EQUAL 2)
        message(FATAL_ERROR "Malformed matrix selection line: ${matrix_line}")
    endif()

    list(GET matrix_parts 0 matrix_id)
    list(GET matrix_parts 1 matrix_is_fast)

    set(verification_test_name "VerificationMatrix_${matrix_id}")
    add_test(
        NAME ${verification_test_name}
        COMMAND ${Python3_EXECUTABLE} ${VERIFY_MATRIX_SCRIPT}
            --matrix-id ${matrix_id}
            --verification-dir ${VERIFICATION_DIR}
            --fracessa-exe $<TARGET_FILE:fracessa>
    )

    if(matrix_is_fast STREQUAL "1")
        set_tests_properties(${verification_test_name} PROPERTIES LABELS "verification;verification_fast")
    else()
        set_tests_properties(${verification_test_name} PROPERTIES LABELS "verification;verification_full_only")
    endif()
endforeach()

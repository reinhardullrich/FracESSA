# Test executables
set(TEST_LIST
    test_matrix
    test_matrix_parser
    test_supports
    test_linear_solver
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

# Register tests with CTest
add_test(NAME MatrixTests COMMAND test_matrix)
add_test(NAME MatrixParserTests COMMAND test_matrix_parser)
add_test(NAME SupportsTests COMMAND test_supports)
add_test(NAME LinearSolverTests COMMAND test_linear_solver)
add_test(NAME CopositivityTests COMMAND test_copositivity)
add_test(NAME Bitset64Tests COMMAND test_bitset64)
add_test(NAME IntegrationTests COMMAND test_integration)
add_test(NAME FractionTests COMMAND test_fraction)

if(TARGET fracessa)
    # CLI parser black-box test
    set(CLI_PARSER_BLACKBOX_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/cli_parser_blackbox.py")
    add_test(
        NAME CliParserBlackBoxTests
        COMMAND ${Python3_EXECUTABLE} ${CLI_PARSER_BLACKBOX_SCRIPT}
            --fracessa-exe $<TARGET_FILE:fracessa>
    )
    set_tests_properties(CliParserBlackBoxTests PROPERTIES LABELS "cli;parser")
endif()

# Test executables
set(TEST_LIST
    test_matrix
    test_linear_solver
    test_bareiss_lu
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
        Threads::Threads
    )
    target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
endforeach()

# Copy test data to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Register tests with CTest
add_test(NAME MatrixTests COMMAND test_matrix)
add_test(NAME LinearSolverTests COMMAND test_linear_solver)
add_test(NAME BareissLUTests COMMAND test_bareiss_lu)
add_test(NAME CopositivityTests COMMAND test_copositivity)
add_test(NAME Bitset64Tests COMMAND test_bitset64)
add_test(NAME IntegrationTests COMMAND test_integration)
add_test(NAME FractionTests COMMAND test_fraction)

# Example of how to compile and link
# target_compile_options
# target_include_directories
# target_link_libraries

add_executable(ZDD_test ZDD_test.cpp)
target_compile_options(ZDD_test PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(ZDD_test PRIVATE ${PROJECT_SOURCE_DIR} ${Trilinos_INCLUDE_DIRS})
target_link_libraries(ZDD_test PRIVATE ${Trilinos_LIBRARIES})
target_link_libraries(ZDD_test PRIVATE ${Trilinos_TPL_LIBRARIES})
target_link_libraries(ZDD_test PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(ZDD_test PRIVATE MPI::MPI_CXX)
add_test(NAME ZDD COMMAND ZDD_test)
# set_tests_properties(ZDD PROPERTIES
#   PASS_REGULAR_EXPRESSION "TestPassed;All ok"
# )

add_executable(TpetraUtil_test TpetraUtil_test.cpp TpetraUtil.cpp)
target_compile_options(TpetraUtil_test PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(TpetraUtil_test PRIVATE ${PROJECT_SOURCE_DIR} ${Trilinos_INCLUDE_DIRS})
target_link_libraries(TpetraUtil_test PRIVATE ${Trilinos_LIBRARIES})
target_link_libraries(TpetraUtil_test PRIVATE ${Trilinos_TPL_LIBRARIES})
target_link_libraries(TpetraUtil_test PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(TpetraUtil_test PRIVATE MPI::MPI_CXX)
add_test(NAME TpetraUtil COMMAND TpetraUtil_test)
# set_tests_properties(TpetraUtil PROPERTIES
#   PASS_REGULAR_EXPRESSION "TestPassed;All ok"
# )

add_executable(Belos_test Belos_test.cpp TpetraUtil.cpp)
target_compile_options(Belos_test PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(Belos_test PRIVATE ${PROJECT_SOURCE_DIR} ${Trilinos_INCLUDE_DIRS})
target_link_libraries(Belos_test PRIVATE ${Trilinos_LIBRARIES})
target_link_libraries(Belos_test PRIVATE ${Trilinos_TPL_LIBRARIES})
target_link_libraries(Belos_test PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(Belos_test PRIVATE MPI::MPI_CXX)
add_test(NAME Belos COMMAND Belos_test)
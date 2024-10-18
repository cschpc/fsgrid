cmake_minimum_required(VERSION 3.18)

find_package(MPI REQUIRED)

set(testname ${PROJECT_NAME}_unit_tests)

add_executable(
    ${testname}
    tools_tests.cpp
    grid_tests.cpp
    coordinates_tests.cpp
)

target_include_directories(${testname} PRIVATE
    "${PROJECT_SOURCE_DIR}/src"
    )

target_link_libraries(
  ${testname}
  GTest::gtest_main
  MPI::MPI_CXX
)

gtest_discover_tests(${testname})

cmake_minimum_required(VERSION 3.18)

set(testname ${PROJECT_NAME}_unit_tests)

add_executable(
    ${testname}
    foo.cpp
)

target_include_directories(${testname} PRIVATE
    "${PROJECT_SOURCE_DIR}"
    )

target_link_libraries(
  ${testname}
  GTest::gtest_main
)

gtest_discover_tests(${testname})

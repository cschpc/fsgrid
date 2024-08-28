cmake_minimum_required(VERSION 3.18)

include(GoogleTest)

set(data_dirname testdata)
set(data_src_dir ${CMAKE_CURRENT_SOURCE_DIR}/${data_dirname})

add_subdirectory(unit_tests)
add_subdirectory(mpi_tests)

add_custom_target(testdata_symlink ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${data_src_dir}" "${CMAKE_CURRENT_BINARY_DIR}/${data_dirname}"
    COMMENT "Add symbolic link from ${data_src_dir} to ${CMAKE_CURRENT_BINARY_DIR}/${data_dirname}"
    VERBATIM
)

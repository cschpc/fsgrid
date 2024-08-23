cmake_minimum_required(VERSION 3.18...3.24)

project(${project_name} VERSION 1.0 DESCRIPTION "Fsgrid testing with ctest and googletest" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Download googletest as a dependency
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        f8d7d77c06936315286eb55f8de22cd23c188571 # release-1.14.0
)
FetchContent_MakeAvailable(googletest)

include(CTest)
add_subdirectory(tests)
enable_testing()

set(gcc "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")

set(debug_gcc "$<AND:$<CONFIG:Debug>,${gcc}>")

set(gcc_warn -Wall;-Werror;-Wextra;-Wconversion;-Wsign-conversion;)
set(gcc_warn ${gcc_warn};-pedantic-errors;-Wcast-qual;-Wwrite-strings;)
set(gcc_warn ${gcc_warn};-Wcast-align=strict;-Wparentheses;)
set(gcc_warn ${gcc_warn};-Wlogical-op;-Wlogical-not-parentheses;)
set(gcc_warn ${gcc_warn};-Wredundant-decls;-Wformat=2;)
set(gcc_warn ${gcc_warn};-Wformat-security;-Wformat-nonliteral;)
set(gcc_warn ${gcc_warn};-Wnull-dereference;-Winit-self;-Wuninitialized;)
set(gcc_warn ${gcc_warn};-Warith-conversion;-Wduplicated-branches;)
set(gcc_warn ${gcc_warn};-Wpointer-arith;-Wundef;)
set(gcc_warn ${gcc_warn};-Wduplicated-cond;-Wformat-signedness;)

set(gcc_deb_opt -Og;)
set(gcc_rel_opt -O3;)

add_compile_options(
    "$<${gcc}:${gcc_warn}>"
    "$<IF:${debug_gcc},${gcc_deb_opt},${gcc_rel_opt}>"
    )

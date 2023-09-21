cmake_minimum_required(VERSION 3.21 FATAL_ERROR)

project(fortran_numeric_lib LANGUAGES Fortran)


# defines targets and sources
add_subdirectory(src)

# contains an "external" library we will link to
# add_subdirectory(external)

# enable testing and define tests
enable_testing()
add_subdirectory(tests)
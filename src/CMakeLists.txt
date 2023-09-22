add_library(fortran_numeric_lib "")

target_sources(fortran_numeric_lib
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/interpolation.f90
    ${CMAKE_CURRENT_LIST_DIR}/list.f90
    )
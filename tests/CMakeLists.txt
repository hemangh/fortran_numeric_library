add_executable(interpolation_test 
                test_interpolation.f90)

target_link_libraries(interpolation_test fortran_numeric_lib)

add_executable(list_test 
                test_list.f90)

target_link_libraries(list_test fortran_numeric_lib)

add_test(
  NAME
    test_interpolation
  COMMAND
    $<TARGET_FILE:interpolation_test>

  )

add_test(
    NAME
      test_list
    COMMAND
      $<TARGET_FILE:list_test>
  
    )
FILE(GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test*.f90)

add_executable(test_aaaaaa test_aaaaaa.f90)
target_link_libraries(test_aaaaaa slatec  )

add_executable(test_arrays_01 test_arrays_01.f90)
target_link_libraries(test_arrays_01 arrays  )

add_executable(test_blas_01 test_blas_01.f90)
target_link_libraries(test_blas_01 blas dqdot)

add_executable(test_blas_02 test_blas_02.f90 chkxer.f90 )
target_link_libraries(test_blas_02 blas )

add_executable(test_blas_03 test_blas_03.f90 chkxer.f90)
target_link_libraries(test_blas_03 blas )

add_executable(test_blas_04 test_blas_04.f90 chkxer.f90)
target_link_libraries(test_blas_04 blas )

add_executable(test_bspline_01 test_bspline_01.f90)
target_link_libraries(test_bspline_01 bspline fitting)

add_executable(test_bspline_02 test_bspline_02.f90)
target_link_libraries(test_bspline_02 bspline fitting)

add_executable(test_dassl_01 test_dassl_01.f90)
target_link_libraries(test_dassl_01 dassl )

add_executable(test_dassl_02 test_dassl_02.f90)
target_link_libraries(test_dassl_02 dassl )

add_executable(test_depac_01 test_depac_01.f90 pass.f90)
target_link_libraries(test_depac_01 depac bvsup)

add_executable(test_depac_02 test_depac_02.f90 pass.f90)
target_link_libraries(test_depac_02 depac bvsup)

add_executable(test_eispack_01 test_eispack_01.f90)
target_link_libraries(test_eispack_01 eispack lineq)

add_executable(test_fishpack_01 test_fishpack_01.f90)
target_link_libraries(test_fishpack_01 fishpack )

add_executable(test_fishpack_02 test_fishpack_02.f90)
target_link_libraries(test_fishpack_02 fishpack )

add_executable(test_fitting_01 test_fitting_01.f90)
target_link_libraries(test_fitting_01 fitting least_squares splp)

add_executable(test_fitting_02 test_fitting_02.f90)
target_link_libraries(test_fitting_02 fitting least_squares splp)

add_executable(test_integration_01 test_integration_01.f90)
target_link_libraries(test_integration_01 integration )

add_executable(test_integration_02 test_integration_02.f90)
target_link_libraries(test_integration_02 integration )

add_executable(test_interpolation_01 test_interpolation_01.f90)
target_link_libraries(test_interpolation_01 interpolation  )

add_executable(test_lineq_01 test_lineq_01.f90)
target_link_libraries(test_lineq_01 lineq )

add_executable(test_lineq_02 test_lineq_02.f90)
target_link_libraries(test_lineq_02 lineq )

add_executable(test_linpack_01 test_linpack_01.f90)
target_link_libraries(test_linpack_01 linpack )

add_executable(test_minpack_01 test_minpack_01.f90 pass.f90)
target_link_libraries(test_minpack_01 minpack blas roots)

add_executable(test_minpack_02 test_minpack_02.f90 pass.f90)
target_link_libraries(test_minpack_02 minpack blas roots)

add_executable(test_minpack_03 test_minpack_03.f90 pass.f90)
target_link_libraries(test_minpack_03 minpack fitting splp least_squares)

add_executable(test_minpack_04 test_minpack_04.f90 pass.f90)
target_link_libraries(test_minpack_04 minpack fitting splp least_squares)

add_executable(test_pchip_01 test_pchip_01.f90 comp.f90)
target_link_libraries(test_pchip_01 pchip bspline special_functions)

add_executable(test_pchip_02 test_pchip_02.f90 comp.f90)
target_link_libraries(test_pchip_02 pchip bspline special_functions)

add_executable(test_quadpack_01 test_quadpack_01.f90)
target_link_libraries(test_quadpack_01 quadpack )

add_executable(test_quadpack_02 test_quadpack_02.f90)
target_link_libraries(test_quadpack_02 quadpack )

add_executable(test_roots_01 test_roots_01.f90 pass.f90)
target_link_libraries(test_roots_01 roots )

add_executable(test_sdrive_01 test_sdrive_01.f90)
target_link_libraries(test_sdrive_01 sdrive )

add_executable(test_sdrive_02 test_sdrive_02.f90)
target_link_libraries(test_sdrive_02 sdrive )

add_executable(test_sdrive_03 test_sdrive_03.f90)
target_link_libraries(test_sdrive_03 sdrive )

add_executable(test_slap_01 test_slap_01.f90 ismpl.f90)
target_link_libraries(test_slap_01 slap special_functions)

add_executable(test_slap_02 test_slap_02.f90 ismpl.f90)
target_link_libraries(test_slap_02 slap special_functions)

add_executable(test_specfun_01 test_specfun_01.f90)
target_link_libraries(test_specfun_01 special_functions )

add_executable(test_specfun_02 test_specfun_02.f90)
target_link_libraries(test_specfun_02 special_functions )

add_executable(test_specfun_03 test_specfun_03.f90)
target_link_libraries(test_specfun_03 special_functions )

add_executable(test_specfun_04 test_specfun_04.f90)
target_link_libraries(test_specfun_04 special_functions integration)

add_executable(test_specfun_05 test_specfun_05.f90)
target_link_libraries(test_specfun_05 special_functions integration)

add_executable(test_specfun_06 test_specfun_06.f90)
target_link_libraries(test_specfun_06 special_functions )

add_executable(test_specfun_07 test_specfun_07.f90)
target_link_libraries(test_specfun_07 special_functions )

add_executable(test_specfun_08 test_specfun_08.f90)
target_link_libraries(test_specfun_08 special_functions )

add_executable(test_specfun_09 test_specfun_09.f90)
target_link_libraries(test_specfun_09 special_functions )

add_executable(test_specfun_10 test_specfun_10.f90)
target_link_libraries(test_specfun_10 special_functions )

add_executable(test_specfun_11 test_specfun_11.f90)
target_link_libraries(test_specfun_11 special_functions )

add_executable(test_specfun_12 test_specfun_12.f90)
target_link_libraries(test_specfun_12 special_functions )

add_executable(test_specfun_13 test_specfun_13.f90)
target_link_libraries(test_specfun_13 special_functions )

add_executable(test_splp_01 test_splp_01.f90 pass.f90)
target_link_libraries(test_splp_01 splp least_squares special_functions)

add_executable(test_splp_02 test_splp_02.f90 pass.f90)
target_link_libraries(test_splp_02 splp least_squares)

add_executable(test_wigner_01 test_wigner_01.f90)
target_link_libraries(test_wigner_01 wigner_coef )

add_executable(test_wigner_02 test_wigner_02.f90)
target_link_libraries(test_wigner_02 wigner_coef )

add_library(get_argument OBJECT get_argument.f90)

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe ${tst})
  target_link_libraries(${exe} $<TARGET_OBJECTS:get_argument>)
  add_test(NAME ${exe} COMMAND ./${exe} 1 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties(${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;FAILED")
endforeach()

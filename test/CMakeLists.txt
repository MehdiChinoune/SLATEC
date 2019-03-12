FILE(GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test*.f90)

add_executable(test_aaaaaa test_aaaaaa.f90)

add_executable(test_specfun_01 test_specfun_01.f90)

add_executable(test_specfun_02 test_specfun_02.f90)

add_executable(test_specfun_03 test_specfun_03.f90)

add_executable(test_specfun_04 test_specfun_04.f90)

add_executable(test_specfun_05 test_specfun_05.f90)

add_executable(test_specfun_06 test_specfun_06.f90)

add_executable(test_specfun_07 test_specfun_07.f90)

add_executable(test_specfun_08 test_specfun_08.f90)

add_executable(test_specfun_09 test_specfun_09.f90)

add_executable(test_specfun_10 test_specfun_10.f90)

add_executable(test_specfun_11 test_specfun_11.f90)

add_executable(test_specfun_12 test_specfun_12.f90)

add_executable(test_specfun_13 test_specfun_13.f90)

add_executable(test_wigner_01 test_wigner_01.f90)

add_executable(test_wigner_02 test_wigner_02.f90)

add_executable(test_blas_01 test_blas_01.f90)

add_executable(test_blas_02 test_blas_02.f90 chkxer.f90 )

add_executable(test_blas_03 test_blas_03.f90 chkxer.f90)

add_executable(test_blas_04 test_blas_04.f90 chkxer.f90)

add_executable(test_lineq_01 test_lineq_01.f90)

add_executable(test_lineq_02 test_lineq_02.f90)

add_executable(test_linpack_01 test_linpack_01.f90)

add_executable(test_eispack_01 test_eispack_01.f90)

add_executable(test_slap_01 test_slap_01.f90 ismpl.f90)

add_executable(test_slap_02 test_slap_02.f90 ismpl.f90)

add_executable(test_fitting_01 test_fitting_01.f90)

add_executable(test_fitting_02 test_fitting_02.f90)

add_executable(test_interpolation_01 test_interpolation_01.f90)

add_executable(test_bspline_01 test_bspline_01.f90)

add_executable(test_bspline_02 test_bspline_02.f90)

add_executable(test_pchip_01 test_pchip_01.f90 comp.f90)

add_executable(test_pchip_02 test_pchip_02.f90 comp.f90)

add_executable(test_roots_01 test_roots_01.f90 pass.f90)

add_executable(test_minpack_01 test_minpack_01.f90 pass.f90)

add_executable(test_minpack_02 test_minpack_02.f90 pass.f90)

add_executable(test_splp_01 test_splp_01.f90 pass.f90)

add_executable(test_splp_02 test_splp_02.f90 pass.f90)

add_executable(test_quadpack_01 test_quadpack_01.f90)

add_executable(test_quadpack_02 test_quadpack_02.f90)

add_executable(test_integration_01 test_integration_01.f90)

add_executable(test_integration_02 test_integration_02.f90)

add_executable(test_depac_01 test_depac_01.f90 pass.f90)

add_executable(test_depac_02 test_depac_02.f90 pass.f90)

add_executable(test_sdrive_01 test_sdrive_01.f90)

add_executable(test_sdrive_02 test_sdrive_02.f90)

add_executable(test_sdrive_03 test_sdrive_03.f90)

add_executable(test_dassl_01 test_dassl_01.f90)

add_executable(test_dassl_02 test_dassl_02.f90)

add_executable(test_fishpack_01 test_fishpack_01.f90)

add_executable(test_fishpack_02 test_fishpack_02.f90)

add_executable(test_minpack_03 test_minpack_03.f90 pass.f90)

add_executable(test_minpack_04 test_minpack_04.f90 pass.f90)

add_executable(test_arrays_01 test_arrays_01.f90)

add_library(get_argument OBJECT get_argument.f90)

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe ${tst})
  target_link_libraries(${exe} slatec get_argument)
  add_test(NAME ${exe} COMMAND ./${exe} 1 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties(${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;FAILED")
endforeach()

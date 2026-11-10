set( tests
  test_aaaaaa.f90
  test_approximation_dp.f90
  test_approximation_minpack_01_dp.f90
  test_approximation_minpack_01_sp.f90
  test_approximation_minpack_02_dp.f90
  test_approximation_minpack_02_sp.f90
  test_approximation_sp.f90
  test_data_handling.f90
  test_diff_integ_dp.f90
  test_diff_integ_eq_dassl_dp.f90
  test_diff_integ_eq_dassl_sp.f90
  test_diff_integ_eq_depac_dp.f90
  test_diff_integ_eq_depac_sp.f90
  test_diff_integ_eq_fishpack.f90
  test_diff_integ_eq_sdrive_dp.f90
  test_diff_integ_eq_sdrive_sp_cx.f90
  test_diff_integ_eq_sdrive_sp.f90
  test_diff_integ_quadpack_dp.f90
  test_diff_integ_quadpack_sp.f90
  test_diff_integ_sp.f90
  test_integ_trans_fftpack.f90
  test_interpolation_bspline_dp.f90
  test_interpolation_bspline_sp.f90
  test_interpolation_pchip_dp.f90
  test_interpolation_pchip_sp.f90
  test_interpolation.f90
  test_linear_dqdot.f90
  test_linear_ge.f90
  test_linear_linpack_cx.f90
  test_linear_ponb.f90
  test_linear_slap_dp.f90
  test_linear_slap_sp.f90
  test_nonlin_eq.f90
  test_optimization_splp_dp.f90
  test_optimization_splp_sp.f90
  test_specfun_bes_dp_cx.f90
  test_specfun_bes_dp.f90
  test_specfun_bes_sp_cx.f90
  test_specfun_bes_sp.f90
  test_specfun_cir_dp.f90
  test_specfun_cir_sp.f90
  test_specfun_fnlib_dp.f90
  test_specfun_fnlib_sp_cx.f90
  test_specfun_fnlib_sp.f90
  test_specfun_psi_dp.f90
  test_specfun_psi_sp.f90
  test_specfun_wigner_dp.f90
  test_specfun_wigner_sp.f90
  test_specfun_xleg_dp.f90
  test_specfun_xleg_sp.f90
)

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/modules )

add_library( common_mod common.f90 )
target_link_libraries( common_mod PRIVATE slatec )

set( common_files get_argument.f90 ismpl.f90 pass.f90 )
set_source_files_properties( ${common_files} PROPERTIES HEADER_FILE_ONLY TRUE )
target_sources( common_mod PRIVATE ${common_files} )

foreach( tst ${tests} )
  STRING( REPLACE ".f90" "" exe ${tst} )
  add_executable( ${exe} ${tst} )
  target_link_libraries( ${exe} common_mod slatec )
  add_test( NAME ${exe} COMMAND ./${exe} 2 WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
  set_tests_properties( ${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Failed;FAILED" )
endforeach()

find_package( PFUNIT )

if( PFUNIT_FOUND )
  add_pfunit_ctest( test_aaaaaa_pf
    TEST_SOURCES test_aaaaaa_pf.pf
    LINK_LIBRARIES slatec
  )

  add_subdirectory( approximation )

  add_subdirectory( data_handling )

  add_subdirectory( diff_integ )
endif()

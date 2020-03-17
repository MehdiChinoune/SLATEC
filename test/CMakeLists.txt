FILE( GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS test*.f90 )

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/modules )

add_library( common_mod SHARED common.f90 )
target_link_libraries( common_mod PRIVATE special_functions )

set( common_files get_argument.f90 ismpl.f90 pass.f90 )
set_source_files_properties( ${common_files} PROPERTIES HEADER_FILE_ONLY TRUE )
target_sources( common_mod PRIVATE ${common_files} )

foreach( tst ${tests} )
  STRING( REPLACE ".f90" "" exe ${tst} )
  add_executable( ${exe} ${tst} )
  target_link_libraries( ${exe} common_mod )
  add_test( NAME ${exe} COMMAND ./${exe} 2 WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
  set_tests_properties( ${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Failed;FAILED" )
endforeach()

target_link_libraries( test_aaaaaa slatec )

target_link_libraries( test_approximation_sp approximation )
target_link_libraries( test_approximation_dp approximation )
target_link_libraries( test_approximation_minpack_01_sp approximation nonlin_eq )
target_link_libraries( test_approximation_minpack_01_dp approximation nonlin_eq )
target_link_libraries( test_approximation_minpack_02_sp approximation interpolation )
target_link_libraries( test_approximation_minpack_02_dp approximation interpolation )

target_link_libraries( test_data_handling data_handling )

target_link_libraries( test_diff_integ_sp diff_integ )
target_link_libraries( test_diff_integ_dp diff_integ )
target_link_libraries( test_diff_integ_quadpack_sp diff_integ )
target_link_libraries( test_diff_integ_quadpack_dp diff_integ )

target_link_libraries( test_diff_integ_eq_dassl_sp diff_integ_eq )
target_link_libraries( test_diff_integ_eq_dassl_dp diff_integ_eq )
target_link_libraries( test_diff_integ_eq_depac_sp diff_integ_eq )
target_link_libraries( test_diff_integ_eq_depac_dp diff_integ_eq )
target_link_libraries( test_diff_integ_eq_fishpack diff_integ_eq )
target_link_libraries( test_diff_integ_eq_sdrive_sp diff_integ_eq )
target_link_libraries( test_diff_integ_eq_sdrive_dp diff_integ_eq )
target_link_libraries( test_diff_integ_eq_sdrive_sp_cx diff_integ_eq )

target_link_libraries( test_integ_trans_fftpack integ_trans )

target_link_libraries( test_interpolation interpolation )
target_link_libraries( test_interpolation_bspline_sp interpolation diff_integ )
target_link_libraries( test_interpolation_bspline_dp interpolation diff_integ )
target_link_libraries( test_interpolation_pchip_sp interpolation special_functions )
target_link_libraries( test_interpolation_pchip_dp interpolation special_functions )

target_link_libraries( test_linear_ge linear )
target_link_libraries( test_linear_ponb linear )
target_link_libraries( test_linear_dqdot linear )
target_link_libraries( test_linear_slap_sp linear special_functions )
target_link_libraries( test_linear_slap_dp linear special_functions )
target_link_libraries( test_linear_linpack_cx linear )

target_link_libraries( test_nonlin_eq nonlin_eq )

target_link_libraries( test_optimization_splp_sp optimization approximation )
target_link_libraries( test_optimization_splp_dp optimization approximation )

target_link_libraries( test_specfun_bes_sp special_functions diff_integ )
target_link_libraries( test_specfun_bes_dp special_functions diff_integ )
target_link_libraries( test_specfun_bes_sp_cx special_functions )
target_link_libraries( test_specfun_bes_dp_cx special_functions )
target_link_libraries( test_specfun_cir_sp special_functions )
target_link_libraries( test_specfun_cir_dp special_functions )
target_link_libraries( test_specfun_fnlib_sp special_functions )
target_link_libraries( test_specfun_fnlib_dp special_functions )
target_link_libraries( test_specfun_fnlib_sp_cx special_functions )
target_link_libraries( test_specfun_psi_sp special_functions )
target_link_libraries( test_specfun_psi_dp special_functions )
target_link_libraries( test_specfun_wigner_sp special_functions )
target_link_libraries( test_specfun_wigner_dp special_functions )
target_link_libraries( test_specfun_xleg_sp special_functions )
target_link_libraries( test_specfun_xleg_dp special_functions )

find_package( PFUNIT )

add_pfunit_test( test_aaaaaa_pf test_aaaaaa_pf.pf "" "" )
target_link_libraries( test_aaaaaa_pf slatec )

add_subdirectory( approximation )

add_subdirectory( data_handling )

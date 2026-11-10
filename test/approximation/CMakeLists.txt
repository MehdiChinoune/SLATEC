set( tests test_dcov_xf.f90
  test_dfc_xf_1.f90
  test_dfc_xf_2.f90
  test_dfc_xf_3.f90
  test_dfc_xf_4.f90
  test_dfc_xf_5.f90
  test_dfc_xf_6.f90
  test_dlsei_xf_1.f90
  test_dlsei_xf_2.f90
  test_dnls1e_xf.f90
  test_dnsqe_xf.f90
  test_dpfitt_xf.f90
  test_fc_xf_1.f90
  test_fc_xf_2.f90
  test_fc_xf_3.f90
  test_fc_xf_4.f90
  test_fc_xf_5.f90
  test_fc_xf_6.f90
  test_lsei_xf_1.f90
  test_lsei_xf_2.f90
  test_pfitqx_xf.f90
  test_scov_xf.f90
  test_snls1e_xf.f90
  test_snsqe_xf.f90
)

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe ${tst})
  add_executable( ${exe} ${tst} )
  target_link_libraries( ${exe} slatec )
  add_test( NAME ${exe} COMMAND ./${exe} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties( ${exe} PROPERTIES WILL_FAIL ON )
endforeach()

set( tests_pf
  test_dfc_pf.pf
  test_dlsei_pf.pf
  test_dnls1q_pf.pf
  test_dnsqe_pf.pf
  test_dpfitt_pf.pf
  test_fc_pf.pf
  test_lsei_pf.pf
  test_pfitqx_pf.pf
  test_snls1q_pf.pf
  test_snsqe_pf.pf
)

foreach( tst ${tests_pf} )
  STRING(REPLACE ".pf" "" exe ${tst})
  add_pfunit_ctest( ${exe}
    TEST_SOURCES ${tst}
    LINK_LIBRARIES slatec
  )
endforeach()

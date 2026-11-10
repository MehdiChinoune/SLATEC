set( tests_pf
  test_diff_integ_dp_pf.pf
  test_diff_integ_quadpack_dp_pf.pf
  test_diff_integ_quadpack_sp_pf.pf
  test_diff_integ_sp_pf.pf
)

foreach( tst ${tests_pf} )
  STRING(REPLACE ".pf" "" exe ${tst})
  add_pfunit_ctest( ${exe}
    TEST_SOURCES ${tst}
    LINK_LIBRARIES slatec
  )
endforeach()

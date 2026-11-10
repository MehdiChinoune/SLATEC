set( tests
  test_dpperm_xf_1.f90
  test_dpperm_xf_2.f90
  test_dpsort_xf_1.f90
  test_dpsort_xf_2.f90
  test_hpperm_xf_1.f90
  test_hpperm_xf_2.f90
  test_hpsort_xf_1.f90
  test_hpsort_xf_2.f90
  test_hpsort_xf_3.f90
  test_hpsort_xf_4.f90
  test_hpsort_xf_5.f90
  test_hpsort_xf_6.f90
  test_ipperm_xf_1.f90
  test_ipperm_xf_2.f90
  test_ipsort_xf_1.f90
  test_ipsort_xf_2.f90
  test_spperm_xf_1.f90
  test_spperm_xf_2.f90
  test_spsort_xf_1.f90
  test_spsort_xf_2.f90
)

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe ${tst})
  add_executable( ${exe} ${tst} )
  target_link_libraries( ${exe} slatec )
  add_test( NAME ${exe} COMMAND ./${exe} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties( ${exe} PROPERTIES WILL_FAIL ON )
endforeach()

set( tests_pf test_data_handling_pf.pf )

foreach( tst ${tests_pf} )
  STRING(REPLACE ".pf" "" exe ${tst})
  add_pfunit_ctest( ${exe}
    TEST_SOURCES ${tst}
    LINK_LIBRARIES slatec
  )
endforeach()

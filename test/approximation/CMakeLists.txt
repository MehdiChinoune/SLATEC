FILE(GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS test*.f90 )

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe ${tst})
  add_executable( ${exe} ${tst} )
  target_link_libraries( ${exe} approximation )
  add_test( NAME ${exe} COMMAND ./${exe} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties( ${exe} PROPERTIES WILL_FAIL ON )
endforeach()

FILE(GLOB tests_pf RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS test*.pf )

foreach( tst ${tests_pf} )
  STRING(REPLACE ".pf" "" exe ${tst})
  add_pfunit_test( ${exe} ${tst} "" "" )
  target_link_libraries( ${exe} approximation )
endforeach()

target_link_libraries( test_dfc_pf interpolation )
target_link_libraries( test_fc_pf interpolation )

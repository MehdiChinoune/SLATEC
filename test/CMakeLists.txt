FILE( GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS test*.f90 )

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

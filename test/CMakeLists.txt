FILE(GLOB tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test*.f90 )

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/modules )

add_library( common_mod SHARED common.f90 )
target_link_libraries( common_mod slatec )

set( common_files get_argument.f90 ismpl.f90 pass.f90 )
set_source_files_properties( ${common_files} PROPERTIES HEADER_FILE_ONLY TRUE )
target_sources( common_mod PRIVATE ${common_files} )

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe ${tst})
  add_executable(${exe} ${tst} )
  target_link_libraries(${exe} common_mod )
  add_test(NAME ${exe} COMMAND ./${exe} 2 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  set_tests_properties(${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Failed;FAILED")
endforeach()

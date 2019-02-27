FILE(GLOB tests test*.f90)

add_executable(test01 test01.f90)

add_executable(test02 test02.f90)

add_executable(test03 test03.f90)

add_executable(test04 test04.f90)

add_executable(test05 test05.f90)

add_executable(test06 test06.f90)

add_executable(test07 test07.f90)

add_executable(test08 test08.f90)

add_executable(test09 test09.f90)

add_executable(test10 test10.f90)

add_executable(test11 test11.f90)

add_executable(test12 test12.f90)

add_executable(test13 test13.f90)

add_executable(test14 test14.f90)

add_executable(test15 test15.f90)

add_executable(test16 test16.f90)

add_executable(test17 test17.f90)

add_executable(test18 test18.f90 chkxer.f90 )

add_executable(test19 test19.f90 chkxer.f90)

add_executable(test20 test20.f90 chkxer.f90)

add_executable(test21 test21.f90)

add_executable(test22 test22.f90)

add_executable(test23 test23.f90)

add_executable(test24 test24.f90)

add_executable(test25 test25.f90 ismpl.f90)

add_executable(test26 test26.f90 ismpl.f90)

add_executable(test27 test27.f90)

add_executable(test28 test28.f90)

add_executable(test29 test29.f90)

add_executable(test30 test30.f90)

add_executable(test31 test31.f90)

add_executable(test32 test32.f90 comp.f90)

add_executable(test33 test33.f90 comp.f90)

add_executable(test34 test34.f90 pass.f90)

add_executable(test35 test35.f90 pass.f90)

add_executable(test36 test36.f90 pass.f90)

add_executable(test37 test37.f90 pass.f90)

add_executable(test38 test38.f90 pass.f90)

add_executable(test39 test39.f90)

add_executable(test40 test40.f90)

add_executable(test41 test41.f90)

add_executable(test42 test42.f90)

add_executable(test43 test43.f90 pass.f90)

add_executable(test44 test44.f90 pass.f90)

add_executable(test45 test45.f90)

add_executable(test46 test46.f90)

add_executable(test47 test47.f90)

add_executable(test48 test48.f90)

add_executable(test49 test49.f90)

add_executable(test50 test50.f90)

add_executable(test51 test51.f90)

add_executable(test52 test52.f90 pass.f90)

add_executable(test53 test53.f90 pass.f90)

add_executable(test54 test54.f90)

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe1 ${tst})
  STRING(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" exe ${exe1})
  target_link_libraries(${exe} slatec)
  if(UNIX)
    add_test(NAME ${exe} COMMAND sh -c "./${exe} < ${CMAKE_CURRENT_SOURCE_DIR}/kprint.txt" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  elseif(WIN32)
    add_test(NAME ${exe} COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} $ENV{ComSpec} /c "${exe} < ${CMAKE_CURRENT_SOURCE_DIR}/kprint.txt")
  endif()
  set_tests_properties(${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;FAILED")
endforeach()

cmake_minimum_required(VERSION 3.7)

project(slatec Fortran)

enable_testing()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -std=f2003 -fimplicit-none" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -fcheck=all -g -Og -std=f2003" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -warn all -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -warn all -g" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 " )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -g" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -g -C" )

endif ()

file(GLOB srcs src/*.f90 fnlib/*.f90 pchip/*.f90 fishfft/*.f90 linear/*.f90 xerror/*.f90 )

add_library(slatec SHARED ${srcs})

add_subdirectory(test)

add_executable(spfuntest spfun.f90 spfunchk.f90)
target_link_libraries(spfuntest slatec)
if (UNIX)
  add_test(NAME spfuntest COMMAND sh -c "./spfuntest <${CMAKE_SOURCE_DIR}/chk/kprint.txt" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
elseif (WIN32)
  add_test(NAME spfuntest COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} $ENV{ComSpec} /c "spfuntest <${CMAKE_SOURCE_DIR}/chk/kprint.txt")
endif ()
set_tests_properties(spfuntest PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;FAILED")

add_executable(slprep slprep.f90 subsid.f90)
target_link_libraries(slprep slatec)

add_executable(sladoc sladoc.f90 subsid.f90)
target_link_libraries(sladoc slatec)

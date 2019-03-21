cmake_minimum_required(VERSION 3.7)

project(slatec Fortran)

set ( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "CodeCoverage" )

enable_testing()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -std=f2003 -fimplicit-none" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Wall -Wno-maybe-uninitialized -fcheck=all -g -Og -std=f2018" )
  set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "--coverage -g -O0" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -warn all -check all -stand f18 -g" )
  set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "-prof-gen:srcpos -prof-dir=${CMAKE_BINARY_DIR}/prof/ -g -O0" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 " )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -g" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -g -C" )

endif ()

add_subdirectory(src)

add_subdirectory(test)

add_subdirectory(doc)

add_executable(spfuntest spfun.f90 spfunchk.f90)
target_link_libraries(spfuntest slatec)

add_test(NAME spfuntest COMMAND ./spfuntest 1 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
set_tests_properties(spfuntest PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;FAILED")

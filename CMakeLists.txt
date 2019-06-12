cmake_minimum_required(VERSION 3.7)

project(slatec Fortran)

set( slatec_MAJOR_VERSION 4 )
set( slatec_MINOR_VERSION 2 )
set( slatec_VERSION  ${slatec_MAJOR_VERSION}.${slatec_MINOR_VERSION} )

set ( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "CodeCoverage" )

enable_testing()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

if(NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE "Release" )
  message(STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}")
endif()

if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Wall -Wimplicit-procedure \
    -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace -g -O0 -std=f2018" )
  set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "--coverage -g -O0" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -warn all -traceback \
    -fpe0 -g -check all -stand f18 -standrad-semantics" )
  set ( CMAKE_Fortran_FLAGS_CODECOVERAGE "-prof-gen:srcpos -prof-dir=${CMAKE_BINARY_DIR}/prof/ -g -O0" )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 " )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -g -Wall " )

elseif ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI" )

  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -Minform=inform -Mallocatable=03 -Mstandard -g -C" )

endif ()

if( "${BLA_VENDOR}" MATCHES "Intel" )
  enable_language( C )
endif()
find_package(BLAS REQUIRED)

find_package(LAPACK REQUIRED)

add_subdirectory(src)

add_subdirectory(test)

cmake_minimum_required(VERSION 3.10)

cmake_policy(SET CMP0069 NEW) # INTERPROCEDURAL_OPTIMIZATION
cmake_policy(SET CMP0074 NEW) # to use <PackageName>_ROOT variables
cmake_policy(SET CMP0083 NEW) # POSITION_INDEPENDENT_CODE

project( slatec
  VERSION 4.2
  DESCRIPTION "A comprehensive software library containing over 1400 general purpose mathematical and statistical routines written in Fortran."
  LANGUAGES Fortran
)

set( slatec_MAJOR_VERSION 4 )
set( slatec_MINOR_VERSION 2 )
set( slatec_VERSION  ${slatec_MAJOR_VERSION}.${slatec_MINOR_VERSION} )

set( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "CodeCoverage" "Fast" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}" )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE "Release" )
  message( STATUS "No CMAKE_BUILD_TYPE selected, default to ${CMAKE_BUILD_TYPE}" )
endif()

if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wimplicit-procedure -Wunused-parameter -std=f2018" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace" )
  set( CMAKE_Fortran_FLAGS_CODECOVERAGE "--coverage -g -O0" )
elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel|IntelLLVM" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -stand f18 -standard-semantics -assume nostd_mod_proc_name" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -fpe-all=0 -traceback" )
  set( CMAKE_Fortran_FLAGS_CODECOVERAGE "-prof-gen:srcpos -prof-dir=${CMAKE_BINARY_DIR}/prof/ -g -O0" )
elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang" )
  set( CMAKE_Fortran_FLAGS_FAST "${CMAKE_Fortran_FLAGS} -Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Wall" )
elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI|NVHPC" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -Mstandard" )
  set( CMAKE_Fortran_FLAGS_FAST "${CMAKE_Fortran_FLAGS} -fast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Minform=inform -C" )
endif()

include( CheckIPOSupported )
check_ipo_supported( RESULT IPO_SUPPORTED OUTPUT IPO_MESSAGE LANGUAGES CXX )
if( IPO_SUPPORTED AND NOT MINGW )
  set( CMAKE_INTERPROCEDURAL_OPTIMIZATION ON )
endif()

include( CheckPIESupported )
check_pie_supported( OUTPUT_VARIABLE PIE_MESSAGE LANGUAGES CXX )
if( CMAKE_CXX_LINK_PIE_SUPPORTED )
  set( CMAKE_POSITION_INDEPENDENT_CODE TRUE )
endif()

if( "${BLA_VENDOR}" MATCHES "Intel" )
  enable_language( C )
endif()
find_package( BLAS REQUIRED )

find_package( LAPACK REQUIRED )

add_subdirectory( src )

enable_testing()
add_subdirectory( test )

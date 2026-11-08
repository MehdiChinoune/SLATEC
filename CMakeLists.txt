cmake_minimum_required(VERSION 3.22)

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
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 14.1 )
    set(_fstd f2023)
  else()
    set(_fstd f2018)
  endif()
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wimplicit-procedure -Wunused-parameter -std=${_fstd}" )
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all -ffpe-trap=zero,overflow,invalid -fbacktrace" )
  set( CMAKE_Fortran_FLAGS_CODECOVERAGE "--coverage -g -O0" )
elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM" )
  if(WIN32)
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /warn:all /stand:f2023 -standard-semantics" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /check:all -fpe-all=0 -traceback" )
  else()
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -warn all -stand f2023 -standard-semantics" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -fpe-all=0 -traceback" )
  endif()
  set( CMAKE_Fortran_FLAGS_FAST "-Ofast" )
  set( CMAKE_Fortran_FLAGS_CODECOVERAGE "-prof-gen:srcpos -prof-dir=${CMAKE_BINARY_DIR}/prof/ -g -O0" )
elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fimplicit-none-ext -std=f2018" )
  set( CMAKE_Fortran_FLAGS_FAST "-O3 -ffast-math -fstack-arrays" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -pedantic -Wall" )
elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mallocatable=03 -Minform=warn -Mstandard" )
  set( CMAKE_Fortran_FLAGS_FAST "${CMAKE_Fortran_FLAGS} -fast" )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Minfo=all -Mneginfo=all -Kieee -Ktrap=fp -C" )
endif()

include( CheckIPOSupported )
check_ipo_supported( RESULT IPO_SUPPORTED OUTPUT IPO_MESSAGE LANGUAGES Fortran )
if( IPO_SUPPORTED AND NOT MINGW )
  set( CMAKE_INTERPROCEDURAL_OPTIMIZATION ON )
endif()

#include( CheckPIESupported )
#check_pie_supported( OUTPUT_VARIABLE PIE_MESSAGE LANGUAGES Fortran )
# objects must be position-independent to go into the shared libslatec; the
# PIE check above only decides whether executables are also linked with -pie
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

if( "${BLA_VENDOR}" MATCHES "Intel" )
  enable_language( C )
endif()
find_package( BLAS REQUIRED )

find_package( LAPACK REQUIRED )

option( BUILD_SHARED_LIBS "Build shared libraries" ON )

include(GNUInstallDirs)

add_subdirectory( src )

enable_testing()
add_subdirectory( test )

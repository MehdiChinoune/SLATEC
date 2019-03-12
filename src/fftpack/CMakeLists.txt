
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

#if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
#  list(REMOVE_ITEM srcs "cfftf.f90" "ezffti.f90" "rffti.f90" "rfftb.f90" )
#endif()

add_library(fftpack SHARED ${srcs})
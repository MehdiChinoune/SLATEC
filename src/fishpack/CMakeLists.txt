
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

#if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
#  list(REMOVE_ITEM srcs "cblktr.f90" )
#endif()

add_library(fishpack SHARED ${srcs})

target_link_libraries(fishpack arrays blas fftpack machine)
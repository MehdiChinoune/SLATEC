
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/ )

file(GLOB mods slatec.f90 linear.f90 xerror.f90 fishfft.f90 fnlib.f90 pchip.f90 src.f90)

add_library(slatec_mod OBJECT ${mods})

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include)
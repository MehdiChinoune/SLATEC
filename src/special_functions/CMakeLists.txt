
add_library( special_functions OBJECT special_functions.f90 )

target_link_libraries( special_functions service )

set_source_files_properties(special_functions.f90 PROPERTIES Fortran_PREPROCESS ON)

include(CheckSourceCompiles)
check_source_compiles( Fortran
  [[
PROGRAM TEST_BES
IMPLICIT NONE
INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,37)
REAL(SP) :: X, Y; Y=BESSEL_J0(X)
END PROGRAM
]]
  HAVE_BESSEL_SP
  SRC_EXT ".f90"
)

if( HAVE_BESSEL_SP )
  set_source_files_properties(special_functions.f90 PROPERTIES COMPILE_DEFINITIONS "HAVE_BESSEL_SP" )
endif()

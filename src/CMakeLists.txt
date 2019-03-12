
#file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS
#  blas/*.f90 bspline/*.f90 bvsup/*.f90 dassl/*.f90 depac/*.f90
#  eispack/*.f90 fftpack/*.f90 fishpack/*.f90 linpack/*.f90 pchip/*.f90
#  quadpack/*.f90 sdrive/*.f90 slap/*.f90 special_functions/*.f90 splp/*.f90
#  machine/*.f90 xerror/*.f90 )

add_library(slatec SHARED "aaaaaa.f90" )


add_subdirectory( blas )
target_link_libraries(slatec blas )

add_subdirectory( bspline )
target_link_libraries(slatec bspline )

add_subdirectory( bvsup )
target_link_libraries(slatec bvsup )

add_subdirectory( dassl )
target_link_libraries(slatec dassl )

add_subdirectory( depac )
target_link_libraries(slatec depac )

add_subdirectory( dqdot )
target_link_libraries(slatec dqdot )

add_subdirectory( eispack )
target_link_libraries(slatec eispack )

add_subdirectory( fftpack )
target_link_libraries(slatec fftpack )

add_subdirectory( fishpack )
target_link_libraries(slatec fishpack )

add_subdirectory( fitting )
target_link_libraries(slatec fitting )

add_subdirectory( integration )
target_link_libraries(slatec integration )

add_subdirectory( interpolation )
target_link_libraries(slatec interpolation )

add_subdirectory( least_squares )
target_link_libraries(slatec least_squares )

add_subdirectory( lineq )
target_link_libraries(slatec lineq )

add_subdirectory( linpack )
target_link_libraries(slatec linpack )

add_subdirectory( machine )
target_link_libraries(slatec machine )

add_subdirectory( minpack )
target_link_libraries(slatec minpack )

add_subdirectory( pchip )
target_link_libraries(slatec pchip )

add_subdirectory( quadpack )
target_link_libraries(slatec quadpack )

add_subdirectory( sdrive )
target_link_libraries(slatec sdrive )

add_subdirectory( roots )
target_link_libraries(slatec roots )

add_subdirectory( slap )
target_link_libraries(slatec slap )

add_subdirectory( special_functions )
target_link_libraries(slatec special_functions )

add_subdirectory( splp )
target_link_libraries(slatec splp )

add_subdirectory( arrays )
target_link_libraries(slatec arrays )

add_subdirectory( xerror )
target_link_libraries(slatec xerror )

add_subdirectory( wigner_coef )
target_link_libraries(slatec wigner_coef )
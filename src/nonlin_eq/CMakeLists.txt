
add_library( nonlin_eq SHARED nonlin_eq.f90 )

target_link_libraries( nonlin_eq linear service ${LAPACK_LIBRARIES} )

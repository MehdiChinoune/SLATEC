
add_library( nonlin_eq OBJECT nonlin_eq.f90 )

target_link_libraries( nonlin_eq linear )


add_library( interpolation OBJECT interpolation.f90 )

target_link_libraries( interpolation service )

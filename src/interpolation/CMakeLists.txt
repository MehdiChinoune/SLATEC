
add_library( interpolation SHARED interpolation.f90 )

target_link_libraries( interpolation service )

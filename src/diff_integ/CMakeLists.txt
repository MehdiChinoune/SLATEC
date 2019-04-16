
add_library( diff_integ SHARED diff_integ.f90 )

target_link_libraries( diff_integ linear interpolation )

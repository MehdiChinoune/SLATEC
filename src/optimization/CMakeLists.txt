
add_library( optimization SHARED optimization.f90 )

target_link_libraries( optimization linear )

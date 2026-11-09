
add_library( optimization OBJECT optimization.f90 )

target_link_libraries( optimization service )

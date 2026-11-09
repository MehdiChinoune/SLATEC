
add_library( approximation OBJECT approximation.f90 )

target_link_libraries( approximation data_handling linear service )

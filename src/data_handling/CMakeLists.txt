
add_library( data_handling SHARED data_handling.f90 )

target_link_libraries( data_handling linear )

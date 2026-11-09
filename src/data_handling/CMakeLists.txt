
add_library( data_handling OBJECT data_handling.f90 )

target_link_libraries( data_handling service )

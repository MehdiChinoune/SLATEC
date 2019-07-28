
add_library( linear SHARED linear.f90 )

target_link_libraries( linear ${LAPACK_LIBRARIES} )

add_dependencies( linear service )

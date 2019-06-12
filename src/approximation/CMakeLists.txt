
add_library( approximation SHARED approximation.f90 )

target_link_libraries( approximation optimization data_handling linear ${BLAS_LIBRARIES} )

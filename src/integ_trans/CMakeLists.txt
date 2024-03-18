
add_library( integ_trans SHARED integ_trans.f90 )

target_link_libraries( integ_trans service )

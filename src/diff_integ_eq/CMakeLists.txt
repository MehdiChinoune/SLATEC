
add_library( diff_integ_eq SHARED diff_integ_eq.f90 )

target_link_libraries( diff_integ_eq data_handling linear integ_trans )


file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 bvsup/*.f90 dassl/*.f90 depac/*.f90 fishpack/*.f90 sdrive/*.f90 )

add_library(diff_integ_eq SHARED ${srcs})

target_link_libraries(diff_integ_eq linear data_handling integ_trans service)
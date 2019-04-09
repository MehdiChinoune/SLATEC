
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 quadpack/*.f90 )

add_library(diff_integ SHARED ${srcs})

target_link_libraries(diff_integ interpolation linear service)
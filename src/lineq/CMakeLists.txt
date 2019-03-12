
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(lineq SHARED ${srcs})

target_link_libraries(lineq blas eispack linpack machine xerror)
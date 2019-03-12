
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(least_squares SHARED ${srcs})

target_link_libraries(least_squares blas splp machine xerror)
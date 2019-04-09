
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(minpack SHARED ${srcs})

target_link_libraries(minpack machine xerror)
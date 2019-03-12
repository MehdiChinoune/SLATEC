
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(ext_exp SHARED ${srcs})

target_link_libraries(ext_exp machine xerror)
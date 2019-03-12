
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(pchip SHARED ${srcs})

target_link_libraries(pchip machine xerror)

file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(slap SHARED ${srcs})

target_link_libraries(slap special_functions blas machine xerror )
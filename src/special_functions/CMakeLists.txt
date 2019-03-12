
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(special_functions SHARED ${srcs})

target_link_libraries(special_functions machine xerror)
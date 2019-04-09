
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS splp/*.f90 )

add_library(optimization SHARED ${srcs})

target_link_libraries(optimization linear service )
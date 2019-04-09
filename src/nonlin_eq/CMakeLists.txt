
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(nonlin_eq SHARED ${srcs})

target_link_libraries(nonlin_eq linear service)

file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 minpack/*.f90 )

add_library(approximation SHARED ${srcs})

target_link_libraries(approximation optimization data_handling linear service )
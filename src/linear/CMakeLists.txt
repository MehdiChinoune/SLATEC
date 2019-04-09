
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 blas/*.f90 dqdot/*.f90 eispack/*.f90 linpack/*.f90 slap/*.f90 )

add_library(linear SHARED ${srcs})

target_link_libraries(linear service )
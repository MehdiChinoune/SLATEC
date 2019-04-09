
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS fftpack/*.f90 )

add_library(integ_trans SHARED ${srcs})
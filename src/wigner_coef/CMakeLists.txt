
file(GLOB srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS *.f90 )

add_library(wigner_coef SHARED ${srcs})

target_link_libraries(wigner_coef machine xerror)

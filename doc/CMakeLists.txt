
add_executable(slprep slprep.f90 subsid.f90)
target_link_libraries(slprep slatec)

add_executable(sladoc sladoc.f90 subsid.f90)
target_link_libraries(sladoc slatec)
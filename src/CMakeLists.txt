
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )

add_library( slatec SHARED slatec.f90 )

file(GLOB_RECURSE srcs CONFIGURE_DEPENDS *.f90 )
list( FILTER srcs EXCLUDE REGEX "slatec.f90" )
set_source_files_properties( ${srcs} PROPERTIES HEADER_FILE_ONLY TRUE )
target_sources( slatec PRIVATE ${srcs} )

add_subdirectory( approximation )
target_link_libraries( slatec approximation )

add_subdirectory( data_handling )
target_link_libraries( slatec data_handling )

add_subdirectory( diff_integ )
target_link_libraries( slatec diff_integ )

add_subdirectory( diff_integ_eq )
target_link_libraries( slatec diff_integ_eq )

add_subdirectory( integ_trans )
target_link_libraries( slatec integ_trans )

add_subdirectory( interpolation )
target_link_libraries( slatec interpolation )

add_subdirectory( linear )
target_link_libraries( slatec linear )

add_subdirectory( nonlin_eq )
target_link_libraries( slatec nonlin_eq )

add_subdirectory( optimization )
target_link_libraries( slatec optimization )

add_subdirectory( service )
target_link_libraries( slatec service )

add_subdirectory( special_functions )
target_link_libraries( slatec special_functions )

install(TARGETS slatec LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules DESTINATION include PATTERN "*_genmod*" EXCLUDE )
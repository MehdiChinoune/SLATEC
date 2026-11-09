
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules )

add_subdirectory( approximation )
add_subdirectory( data_handling )
add_subdirectory( diff_integ )
add_subdirectory( diff_integ_eq )
add_subdirectory( integ_trans )
add_subdirectory( interpolation )
add_subdirectory( linear )
add_subdirectory( nonlin_eq )
add_subdirectory( optimization )
add_subdirectory( service )
add_subdirectory( special_functions )

add_library( slatec
  aaaaaa.f90
  slatec.f90
  $<TARGET_OBJECTS:approximation>
  $<TARGET_OBJECTS:data_handling>
  $<TARGET_OBJECTS:diff_integ>
  $<TARGET_OBJECTS:diff_integ_eq>
  $<TARGET_OBJECTS:integ_trans>
  $<TARGET_OBJECTS:interpolation>
  $<TARGET_OBJECTS:linear>
  $<TARGET_OBJECTS:nonlin_eq>
  $<TARGET_OBJECTS:optimization>
  $<TARGET_OBJECTS:service>
  $<TARGET_OBJECTS:special_functions>
)

target_link_libraries(slatec PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_source_files_properties( aaaaaa.f90 PROPERTIES HEADER_FILE_ONLY TRUE )
target_sources( slatec PRIVATE aaaaaa.f90 )

install( TARGETS slatec
         RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install( DIRECTORY ${CMAKE_BINARY_DIR}/modules/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

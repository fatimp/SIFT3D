add_library (sift3D SHARED sift.c imutil.c nifti.c)
target_link_libraries (sift3D PRIVATE m z ${LAPACK_LIBRARIES})
target_include_directories(sift3D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if (NIFTI_FOUND)
  target_compile_definitions(sift3D PRIVATE "SIFT3D_WITH_NIFTI")
  target_include_directories(sift3D PRIVATE ${NIFTI_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIR})
  target_link_libraries (sift3D PRIVATE ${NIFTI_LIBRARIES} ${ZLIB_LIBRARIES})
endif (NIFTI_FOUND)

set_target_properties (sift3D PROPERTIES VERSION ${SIFT3D_VERSION}
  SOVERSION ${SIFT3D_VERSION_MAJOR}
  C_VISIBILITY_PRESET hidden
)

install (FILES imtypes.h imutil.h sift.h DESTINATION include/sift3d)
install (TARGETS sift3D LIBRARY DESTINATION lib)

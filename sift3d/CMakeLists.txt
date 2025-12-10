add_library (sift3D SHARED sift.c imutil.c nifti.c)
target_link_libraries (sift3D PRIVATE m ${LAPACK_LIBRARIES})

if (NIFTI_FOUND)
  target_compile_definitions(sift3D PRIVATE "SIFT3D_WITH_NIFTI")
  target_include_directories(sift3D PRIVATE ${NIFTI_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIR})
  target_link_libraries (sift3D PRIVATE ${NIFTI_LIBRARIES} ${ZLIB_LIBRARIES})
endif (NIFTI_FOUND)

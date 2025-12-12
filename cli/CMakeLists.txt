add_executable (kpSift3D kpSift3D.c)

target_include_directories (kpSift3D PRIVATE ${CMAKE_SOURCE_DIR}/sift3d)
target_link_libraries  (kpSift3D sift3D)

install (TARGETS kpSift3D RUNTIME DESTINATION bin)

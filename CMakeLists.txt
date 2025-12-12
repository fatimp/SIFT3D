cmake_minimum_required (VERSION 3.10)

project (SIFT3D)

set (SIFT3D_VERSION_MAJOR 2)
set (SIFT3D_VERSION_MINOR 0)
set (SIFT3D_VERSION ${SIFT3D_VERSION_MAJOR}.${SIFT3D_VERSION_MINOR})

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package (OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

find_package (LAPACK REQUIRED)
find_package (NIFTI)

if (NIFTI_FOUND)
  find_package (ZLIB REQUIRED)
endif(NIFTI_FOUND)

add_subdirectory (sift3d)
add_subdirectory (cli)

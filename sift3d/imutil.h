/**
 * @file imutil.h
 * @brief Manupulations with images.
 */

#ifndef __SIFT3D_IMUTIL_H__
#define __SIFT3D_IMUTIL_H__

#include <imtypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Returned by sift3d_read_image() when this file type is not
 * supported.
 */
#define SIFT3D_UNSUPPORTED_FILE_TYPE 2

/**
 * Returned by sift3d_read_image() if the library was compiled without
 * a dependency needed to read this file (e.g. nifticlib).
 */
#define SIFT3D_WRAPPER_NOT_COMPILED  3

// Images

/**
 * @brief Make an image with specified dimensions.
 *
 * The image must be freed with sift3d_free_image() when not needed.
 *
 * @param nx Dimension 0
 * @param ny Dimension 1
 * @param nz Dimension 2
 * @param nc Number of channels. Must be 1.
 */
SIFT3D_EXPORT sift3d_image*
sift3d_make_image (const int nx, const int ny, const int nz, const int nc);

/**
 * @brief Destroy an image.
 */
SIFT3D_EXPORT void
sift3d_free_image (sift3d_image*);

/**
 * @brief Read an image from a file.
 *
 * The image must be freed with sift3d_free_image() when not needed.
 * Currently only NifTi-1 images are supported (.nii or .nii.gz).
 */
SIFT3D_EXPORT sift3d_image*
sift3d_read_image (const char *path);

/**
 * @brief Get a pointer to the image's data
 *
 * The user is supposed to fill the image using this pointer after the
 * image was created with sift3d_make_image(). **NB: Due to a strange
 * whim of the author the image is in column-major order!**
 */
SIFT3D_EXPORT float*
sift3d_image_data (const sift3d_image*);

// Matrices

/**
 * @brief Create a matrix.
 *
 * SIFT3D returns descriptors in arrays wrapped in this type. The
 * actual array allocations are performed under the hood. Users of
 * this library must free the matrix and its underlying array with
 * sift3d_free_mat_rm when it's not needed anymore.
 */
SIFT3D_EXPORT sift3d_mat_rm*
sift3d_make_mat_rm();

/**
 * @brief Destroy a matrix.
 */
SIFT3D_EXPORT void
sift3d_free_mat_rm (sift3d_mat_rm*);

/**
 * @brief Get a pointer to the data array of a matrix.
 *
 * An actual type of the data can be queried with
 * sift3d_mat_rm_type().
 */
SIFT3D_EXPORT void*
sift3d_mat_rm_data (sift3d_mat_rm*);

/**
 * @brief Get dimensions of a matrix.
 *
 * @param num_cols Points to an area where the number of columns is
 * written (if not `NULL`).
 * @param num_rows Points to an area where the number of rows is
 * written (if not `NULL`).
 */
SIFT3D_EXPORT void
sift3d_mat_rm_dimensions (const sift3d_mat_rm*, int *num_cols, int *num_rows);

/**
 * @brief Get the type of elements in a matrix.
 */
SIFT3D_EXPORT sift3d_mat_type
sift3d_mat_rm_type (const sift3d_mat_rm*);

#ifdef __cplusplus
}
#endif

#endif

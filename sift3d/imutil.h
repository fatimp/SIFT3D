#ifndef __SIFT3D_IMUTIL_H__
#define __SIFT3D_IMUTIL_H__

#include <imtypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Extra return codes for this module */
#define SIFT3D_UNSUPPORTED_FILE_TYPE 2 /* The file type is not supported */
#define SIFT3D_WRAPPER_NOT_COMPILED  3 /* Whatever */

// Images
SIFT3D_EXPORT sift3d_image*
sift3d_make_image (const int nx, const int ny, const int nz, const int nc);

SIFT3D_EXPORT void
sift3d_free_image (sift3d_image*);

SIFT3D_EXPORT sift3d_image*
sift3d_read_image (const char *path);

SIFT3D_EXPORT float*
sift3d_image_data (const sift3d_image*);

// Matrices
SIFT3D_EXPORT sift3d_mat_rm*
sift3d_make_mat_rm();

SIFT3D_EXPORT void
sift3d_free_mat_rm (sift3d_mat_rm*);

SIFT3D_EXPORT void*
sift3d_mat_rm_data (sift3d_mat_rm*);

SIFT3D_EXPORT void
sift3d_mat_rm_dimensions (const sift3d_mat_rm*, int *num_cols, int *num_rows);

SIFT3D_EXPORT sift3d_mat_type
sift3d_mat_rm_type (const sift3d_mat_rm*);

#ifdef __cplusplus
}
#endif

#endif

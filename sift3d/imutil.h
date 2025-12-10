/* -----------------------------------------------------------------------------
 * imutil.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header for imutil.c
 * -----------------------------------------------------------------------------
 */

#include "imtypes.h"

#ifndef _IMUTIL_H
#define _IMUTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Extra return codes for this module */
#define SIFT3D_UNSUPPORTED_FILE_TYPE 2 /* The file type is not supported */
#define SIFT3D_WRAPPER_NOT_COMPILED  3 /* Whatever */

/* Externally-visible routines */
void *SIFT3D_safe_realloc(void *ptr, size_t size);

void init_Mesh(sift3d_mesh * const mesh);

void cleanup_Mesh(sift3d_mesh * const mesh);

int init_Mat_rm(sift3d_mat_rm *const mat, const int num_rows, const int num_cols,
                const sift3d_mat_rm_type type, const int set_zero);

int init_Mat_rm_p(sift3d_mat_rm *const mat, const void *const p, const int num_rows, 
                  const int num_cols, const sift3d_mat_rm_type type, 
                  const int set_zero);

int set_Mat_rm_zero(sift3d_mat_rm *mat);

int copy_Mat_rm(const sift3d_mat_rm *const src, sift3d_mat_rm *const dst);

int resize_Mat_rm(sift3d_mat_rm *const mat); 

int eigen_Mat_rm(sift3d_mat_rm *A, sift3d_mat_rm *Q, sift3d_mat_rm *L);

int transpose_Mat_rm(const sift3d_mat_rm *const src, sift3d_mat_rm *const dst);

int zero_Mat_rm(sift3d_mat_rm *const mat);
 
int identity_Mat_rm(const int n, sift3d_mat_rm *const mat);

void cleanup_Mat_rm(sift3d_mat_rm *mat);

sift3d_im_format im_get_format(const char *path);

int im_read(const char *path, sift3d_image *const im);

int im_write(const char *path, const sift3d_image *const im);

int write_Mat_rm(const char *path, const sift3d_mat_rm *const mat);

int init_im_with_dims(sift3d_image *const im, const int nx, const int ny, const int nz,
                      const int nc);

int im_copy_dims(const sift3d_image *const src, sift3d_image *dst);

int im_copy_data(const sift3d_image *const src, sift3d_image *const dst);

void im_free(sift3d_image *im);

int im_downsample_2x(const sift3d_image *const src, sift3d_image *const dst);

int im_permute(const sift3d_image *const src, const int dim1, const int dim2, 
               sift3d_image *const dst);

void im_default_stride(sift3d_image *const im);

int im_resize(sift3d_image *const im);

float im_max_abs(const sift3d_image *const im);

void im_scale(const sift3d_image *const im);

int im_subtract(sift3d_image *src1, sift3d_image *src2, sift3d_image *dst);

void im_zero(sift3d_image *im);

void init_im(sift3d_image *const im);

int init_Gauss_filter(sift3d_gauss_filter *const gauss, const double sigma,
                      const int dim);

int init_Gauss_incremental_filter(sift3d_gauss_filter *const gauss,
                                  const double s_cur, const double s_next,
                                  const int dim);

int init_Sep_FIR_filter(sift3d_sep_fir_filter *const f, const int dim, const int width,
                        const float *const kernel, const int symmetric);

int apply_Sep_FIR_filter(const sift3d_image *const src, sift3d_image *const dst,
                         sift3d_sep_fir_filter *const f, const double unit);

void cleanup_Sep_FIR_filter(sift3d_sep_fir_filter *const f);

void cleanup_Gauss_filter(sift3d_gauss_filter *gauss);

void init_GSS_filters(sift3d_gss_filters *const gss);

int make_gss(sift3d_gss_filters *const gss, const sift3d_pyramid *const pyr);

void cleanup_GSS_filters(sift3d_gss_filters *const gss);

void init_Pyramid(sift3d_pyramid *const pyr);

int resize_Pyramid(const sift3d_image *const im, const int first_level, 
                   const unsigned int num_kp_levels, const unsigned int num_levels,
                   const int first_octave, const unsigned int num_octaves, 
                   sift3d_pyramid *const pyr);

int set_scales_Pyramid(const double sigma0, const double sigma_n, 
                       sift3d_pyramid *const pyr);

void cleanup_Pyramid(sift3d_pyramid *const pyr);

void init_Slab(sift3d_slab *const slab);

void cleanup_Slab(sift3d_slab *const slab);

int resize_Slab(sift3d_slab *slab, int num, size_t size);

#ifdef __cplusplus
}
#endif

#endif

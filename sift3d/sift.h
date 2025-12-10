/* -----------------------------------------------------------------------------
 * sift.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header for sift.c
 * -----------------------------------------------------------------------------
 */

#include "imtypes.h"

#ifndef _SIFT_H
#define _SIFT_H

#ifdef __cplusplus
extern "C" {
#endif

void init_Keypoint_store(sift3d_keypoint_store *const kp);

int init_Keypoint(sift3d_keypoint *const key);

int resize_Keypoint_store(sift3d_keypoint_store *const kp, const size_t num);

int copy_Keypoint(const sift3d_keypoint *const src, sift3d_keypoint *const dst);

void cleanup_Keypoint_store(sift3d_keypoint_store *const kp);

void init_SIFT3D_Descriptor_store(sift3d_descriptor_store *const desc);

void cleanup_SIFT3D_Descriptor_store(sift3d_descriptor_store *const desc);

int set_peak_thresh_SIFT3D(sift3d_detector *const sift3d,
                           const double peak_thresh);

int set_corner_thresh_SIFT3D(sift3d_detector *const sift3d,
                             const double corner_thresh);

int set_num_kp_levels_SIFT3D(sift3d_detector *const sift3d,
                             const unsigned int num_kp_levels);

int set_sigma_n_SIFT3D(sift3d_detector *const sift3d,
                       const double sigma_n);

int set_sigma0_SIFT3D(sift3d_detector *const sift3d,
                      const double sigma_n);

int init_SIFT3D(sift3d_detector *sift3d);

void cleanup_SIFT3D(sift3d_detector *const sift3d);

int SIFT3D_detect_keypoints(sift3d_detector *const sift3d, const sift3d_image *const im,
                            sift3d_keypoint_store *const kp);

int SIFT3D_have_gpyr(const sift3d_detector *const sift3d);

int SIFT3D_extract_descriptors(sift3d_detector *const sift3d, 
                               const sift3d_keypoint_store *const kp, 
                               sift3d_descriptor_store *const desc);

int Keypoint_store_to_Mat_rm(const sift3d_keypoint_store *const kp, sift3d_mat_rm *const mat);

int SIFT3D_Descriptor_coords_to_Mat_rm(
    const sift3d_descriptor_store *const store, 
    sift3d_mat_rm *const mat);

int SIFT3D_Descriptor_store_to_Mat_rm(const sift3d_descriptor_store *const store, 
                                      sift3d_mat_rm *const mat);

int Mat_rm_to_SIFT3D_Descriptor_store(const sift3d_mat_rm *const mat, 
                                      sift3d_descriptor_store *const store);

int write_Keypoint_store(const char *path, const sift3d_keypoint_store *const kp);

int write_SIFT3D_Descriptor_store(const char *path, 
                                  const sift3d_descriptor_store *const desc);

#ifdef __cplusplus
}
#endif

#endif

#ifndef __SIFT3D__H__
#define __SIFT3D__H__

#include <imtypes.h>

#ifdef __cplusplus
extern "C" {
#endif

// Detector

SIFT3D_EXPORT sift3d_detector*
sift3d_make_detector();

SIFT3D_EXPORT void
sift3d_free_detector(sift3d_detector*);

SIFT3D_EXPORT int
sift3d_detector_set_peak_thresh(sift3d_detector *const, const double);

SIFT3D_EXPORT int
sift3d_detector_set_corner_thresh(sift3d_detector *const, const double);

SIFT3D_EXPORT int
sift3d_detector_set_num_kp_levels(sift3d_detector *const, const unsigned int);

SIFT3D_EXPORT int
sift3d_detector_set_sigma_n(sift3d_detector *const, const double);

SIFT3D_EXPORT int
sift3d_detector_set_sigma0(sift3d_detector *const, const double);

SIFT3D_EXPORT int
sift3d_detect_keypoints(sift3d_detector *const,
                        const sift3d_image *const,
                        sift3d_keypoint_store *const);

SIFT3D_EXPORT int
sift3d_detector_has_gpyr(const sift3d_detector *const);

SIFT3D_EXPORT int
sift3d_extract_descriptors(sift3d_detector *const,
                           const sift3d_keypoint_store *const,
                           sift3d_descriptor_store *const);

// Keypoint store

SIFT3D_EXPORT sift3d_keypoint_store*
sift3d_make_keypoint_store();

SIFT3D_EXPORT void
sift3d_free_keypoint_store(sift3d_keypoint_store*);

SIFT3D_EXPORT int
sift3d_keypoint_store_to_mat_rm(const sift3d_keypoint_store *const,
                                sift3d_mat_rm *const);

SIFT3D_EXPORT int
sift3d_keypoint_store_save(const char *path,
                           const sift3d_keypoint_store *const);

// Descriptor store

SIFT3D_EXPORT sift3d_descriptor_store*
sift3d_make_descriptor_store();

SIFT3D_EXPORT void
sift3d_free_descriptor_store(sift3d_descriptor_store*);

SIFT3D_EXPORT int
sift3d_descriptor_store_save(const char *path,
                             const sift3d_descriptor_store *const);

SIFT3D_EXPORT int
sift3d_descriptor_store_to_mat_rm(const sift3d_descriptor_store *const,
                                  sift3d_mat_rm *const);

#ifdef __cplusplus
}
#endif

#endif

/**
 * @file sift.h
 * @brief Detection of keypoints and calculation of descriptors.
 */

#ifndef __SIFT3D__H__
#define __SIFT3D__H__

#include <imtypes.h>

#ifdef __cplusplus
extern "C" {
#endif

// Detector

/**
 * @brief Create a detector.
 *
 * Detector is an object responsible for detection of keypoints and
 * calculation of descriptors. It must be freed with
 * sift3d_free_detector() when it's not needed anymore.
 */
SIFT3D_EXPORT sift3d_detector*
sift3d_make_detector();

/**
 * @brief Destroy a detector.
 */
SIFT3D_EXPORT void
sift3d_free_detector(sift3d_detector*);

/**
 * @brief Set the relative DoG peak threshold.
 *
 * Keypoints which are weaker than this threshold are discarded.
 * The interval for this parameter: [0, 1].
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_detector_set_peak_thresh(sift3d_detector *const, const double);

/**
 * @brief Set the corner threshold.
 *
 * TODO: What is the corner threshold?
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_detector_set_corner_thresh(sift3d_detector *const, const double);

/**
 * @brief Set the number of images in a DoG pyramid per octave.
 *
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_detector_set_num_kp_levels(sift3d_detector *const, const unsigned int);

/**
 * @brief Set the nominal scale parameter of the input data.
 *
 * The parameter must be non-negative.
 *
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_detector_set_sigma_n(sift3d_detector *const, const double);

/**
 * @brief Sets the scale parameter of the first level of octave 0.
 *
 * The parameter must be non-negative.
 *
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_detector_set_sigma0(sift3d_detector *const, const double);

/**
 * @brief Detect keypoints in the image.
 *
 * This function detects coordinates of keypoints and their
 * orientation. The image is not needed anymore after this function
 * returns.
 *
 * @param detector A freshly created detector
 * @param image An image to work with
 * @param store A freshly created keypoint store
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_detect_keypoints(sift3d_detector       *const detector,
                        const sift3d_image    *const image,
                        sift3d_keypoint_store *const store);

/**
 * @brief Extract descriptors from the keypoints.
 *
 * This function extract descriptors from keypoints previously
 * detected with sift3d_detect_keypoints().
 *
 * @param detector The detector used to extract keypoints
 * @param kp_store The keypoint store with keypoints
 * @param desc_store A freshly created descriptor store
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_extract_descriptors(sift3d_detector             *const detector,
                           const sift3d_keypoint_store *const kp_store,
                           sift3d_descriptor_store     *const desc_store);

// Keypoint store

/**
 * @brief Create an empty keypoint store.
 *
 * The created store must be destroyed with
 * sift3d_free_keypoint_store() when it's not needed.
 */
SIFT3D_EXPORT sift3d_keypoint_store*
sift3d_make_keypoint_store();

/**
 * @brief Destroy a keypoint store.
 */
SIFT3D_EXPORT void
sift3d_free_keypoint_store(sift3d_keypoint_store*);

/**
 * @brief Copy keypoints to a matrix.
 *
 * The matrix is resized to dimensions Nx3 where N is the number of
 * keypoints. Each row of the matrix contains coordinates of a
 * keypoint.
 *
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_keypoint_store_to_mat_rm(const sift3d_keypoint_store *const,
                                sift3d_mat_rm *const);

/**
 * @brief Save keypoints to a file.
 *
 * Keypoints are saved in CSV format. Uncompressed and gzipped output
 * is supported.
 *
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_keypoint_store_save(const char *path,
                           const sift3d_keypoint_store *const);

/**
 * @brief Sort keypoints by strength, possibly removing weak keypoints.
 *
 * The keypoints are sorted in descending order, the strongest (the
 * most stable) keypoint appearing first.
 *
 * @param limit If not zero, only at most this amount of the strongest
 * keypoints are kept in the store.
 */
SIFT3D_EXPORT void
sift3d_keypoint_store_sort_by_strength (sift3d_keypoint_store *const, int limit);

// Descriptor store

/**
 * @brief Create an empty descriptor store.
 *
 * The created store must be destroyed with
 * sift3d_free_descriptor_store() when it's not needed.
 */
SIFT3D_EXPORT sift3d_descriptor_store*
sift3d_make_descriptor_store();

/**
 * @brief Destroy a descriptor store.
 */
SIFT3D_EXPORT void
sift3d_free_descriptor_store(sift3d_descriptor_store*);

/**
 * @brief Save descriptors to a file.
 *
 * Descriptors are saved in CSV format. Uncompressed and gzipped
 * output is supported.
 *
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_descriptor_store_save(const char *path,
                             const sift3d_descriptor_store *const);

/**
 * @brief Copy descriptors to a matrix.
 *
 * The matrix is resized to dimensions Nx771 where N is the number of
 * keypoints. Each row of the matrix contains coordinates of a
 * keypoint (the first 3 elements) + 768 elements of a descriptor (the
 * rest of a row).
 *
 * @return `SIFT_SUCCESS` on success, `SIFT_FAILURE` on failure
 */
SIFT3D_EXPORT int
sift3d_descriptor_store_to_mat_rm(const sift3d_descriptor_store *const,
                                  sift3d_mat_rm *const);

#ifdef __cplusplus
}
#endif

#endif

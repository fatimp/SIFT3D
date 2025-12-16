/**
 * @file imtypes.h
 * @brief Types and constants
*/

#ifndef __SIFT3D_IMTYPES_H__
#define __SIFT3D_IMTYPES_H__

#ifdef __cplusplus
extern "C" {
#endif

// Not exactly belongs here, but I don't know a better place
#define SIFT3D_EXPORT __attribute__((visibility ("default")))

// Return codes

/**
 * @brief This core is returned on success by SIFT3D functions.
 */
#define SIFT3D_SUCCESS 0
/**
 * @brief This core is returned on failure by SIFT3D functions.
 */
#define SIFT3D_FAILURE -1

// Truth values
#define SIFT3D_TRUE 1
#define SIFT3D_FALSE 0

typedef struct _sift3d_detector sift3d_detector;
typedef struct _sift3d_keypoint_store sift3d_keypoint_store;
typedef struct _sift3d_descriptor_store sift3d_descriptor_store;
typedef struct _sift3d_image sift3d_image;
typedef struct _sift3d_mat_rm sift3d_mat_rm;

/**
 * Possible data types for matrix elements
 */
typedef enum {
    SIFT3D_DOUBLE,
    SIFT3D_FLOAT,
    SIFT3D_INT
} sift3d_mat_type;

#ifdef __cplusplus
}
#endif

#endif

/* -----------------------------------------------------------------------------
 * imtypes.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This header contains data type definitions.
 * -----------------------------------------------------------------------------
 */

#include <stdlib.h>

#ifndef _IMTYPES_H
#define _IMTYPES_H

#ifdef __cplusplus
extern "C" {
#endif

// Return codes
#define SIFT3D_SUCCESS 0
#define SIFT3D_FAILURE -1

// Truth values
#define SIFT3D_TRUE 1
#define SIFT3D_FALSE 0

// Platform types
#if _Win16 == 1 || _WIN32 == 1 || _WIN64 == 1 ||    \
    defined __WIN32__ || defined __TOS_WIN__ ||     \
    defined __WINDOWS__ && !defined _WINDOWS
#define _WINDOWS
#endif
#if (defined(__MINGW32__) || defined(__MINGW64__)) && defined(_WINDOWS) && \
    !defined _MINGW_WINDOWS
#define _MINGW_WINDOWS
#endif

/* File separator character */
#ifdef _WINDOWS
#define SIFT3D_FILE_SEP '\\'
#else
#define SIFT3D_FILE_SEP '/'
#endif

/* Parameters */
#define NBINS_AZ 8      // Number of bins for azimuthal angles
#define NBINS_PO 4      // Number of bins for polar angles
#define NHIST_PER_DIM 4 // Number of SIFT descriptor histograms per dimension 
#define ICOS_HIST           // Icosahedral gradient histogram

/* Constants */
#define IM_NDIMS 3 // Number of dimensions in an Image
#define ICOS_NFACES 20 // Number of faces in an icosahedron
#define ICOS_NVERT 12 // Number of vertices in an icosahedron

/* Derived constants */
#define DESC_NUM_TOTAL_HIST (NHIST_PER_DIM * NHIST_PER_DIM * NHIST_PER_DIM)
#define DESC_NUMEL (DESC_NUM_TOTAL_HIST * HIST_NUMEL)

// The number of elements in a gradient histogram
#ifdef ICOS_HIST
#define HIST_NUMEL (ICOS_NVERT)
#else
#define HIST_NUMEL (NBINS_AZ * NBINS_PO)
#endif

/* Supported image file formats */
typedef enum {
    ANALYZE, /* Analyze */
    NIFTI, /* NIFTI-1 */ 
    UNKNOWN, /* Not one of the known extensions */
    FILE_ERROR /* Error occurred in determining the format */
} sift3d_im_format;

/* Possible data types for matrix elements */ 
typedef enum {
    SIFT3D_DOUBLE,
    SIFT3D_FLOAT,
    SIFT3D_INT
} sift3d_mat_rm_type;

/* Struct to hold a dense matrix in row-major order */
typedef struct {
    union {
        double *data_double;
        float  *data_float;
        int    *data_int;
    } u;

    size_t size;        // Size of the buffer, in bytes
    int num_cols;           // Number of columns
    int num_rows;           // Number of rows
    int static_mem;         // Flag for statically-allocated memory
    sift3d_mat_rm_type type;       // DOUBLE, FLOAT, or INT
} sift3d_mat_rm;

/* Struct to hold image data. The image is a rectangular prism, 
 * where the bottom-left corner is [0 0 0], the x-stride is 1,
 * the y-stride is the width in x, and the z-stride is the
 * size of an xy plane. For convenience use the macros IM_GET_IDX, 
 * IM_GET_VOX, and IM_SET_VOX to manipulate this struct. */
typedef struct {
    float *data;        // Raster of voxel values ~16MB
    double s;       // scale-space location
    size_t size;        // Total size in pixels
    int nx, ny, nz;     // Dimensions in x, y, and z
    double ux, uy, uz;  // Real world dimensions in x, y, and z
    size_t xs, ys, zs;      // Stride in x, y, and z
    int nc;                 // The number of channels
} sift3d_image;

/* Holds separable FIR filters and programs to apply them */
typedef struct {
    float *kernel;  // filter weights
    int dim;    // dimensionality, e.g. 3 for MRI
    int width;  // number of weights
    int symmetric;  // enable symmetric optimizations: FALSE or TRUE
} sift3d_sep_fir_filter;

/* Holds Gaussian filters */
typedef struct  {
    double sigma;
    sift3d_sep_fir_filter f;
} sift3d_gauss_filter;

/* Holds Gaussian Scale-Space filters */
typedef struct {
    sift3d_gauss_filter first_gauss;   // Used on the very first blur
    sift3d_gauss_filter *gauss_octave; // Array of kernels for one octave
    int num_filters;        // Number of filters for one octave
    int first_level;                // Index of the first scale level
} sift3d_gss_filters;

/* Struct to hold a scale-space image pyramid */
typedef struct {
    // Levels in all octaves
    sift3d_image *levels;

    // Scale-space parameters
    double sigma_n;
    double sigma0;
    int num_kp_levels;

    // Indexing information -- see immacros.h
    int first_octave;
    int num_octaves;
    int first_level;
    int num_levels;
} sift3d_pyramid;

/* Struct defining a vector in spherical coordinates */
typedef struct {
    float mag;  // Magnitude
    float po;   // Polar angle, [0, pi)
    float az;   // Azimuth angle, [0, 2pi)
} sift3d_svec;

/* Struct defining a vector in Cartesian coordinates */
typedef struct {
    float x;
    float y;
    float z;
} sift3d_cvec;

/* Slab allocation struct */
typedef struct {
    void *buf;          // Buffer
    size_t num;         // Number of elements currently in buffer
    size_t buf_size;    // Buffer capactiy, in bytes
} sift3d_slab;

/* Struct defining a keypoint in 3D space. */
typedef struct {
    float r_data[IM_NDIMS * IM_NDIMS];  // Memory for matrix R, do not use this
    sift3d_mat_rm R;               // Rotation matrix into Keypoint space
    double xd, yd, zd;          // sub-pixel x, y, z
    double  sd;             // absolute scale
    int o, s;                   // pyramid indices
} sift3d_keypoint;

/* Struct to hold keypoints */
typedef struct {
    sift3d_keypoint *buf;
    sift3d_slab slab;
    int nx, ny, nz;     // dimensions of first octave
} sift3d_keypoint_store;

/* Struct defining an orientation histogram in
 * spherical coordinates. */
typedef struct {
    float bins[HIST_NUMEL];
} sift3d_hist;

/* Triangle */
typedef struct {
    sift3d_cvec v[3]; // Vertices
    int idx[3]; // Index of each vertex in the solid
} sift3d_tri;

/* Triangle mesh */
typedef struct {
    sift3d_tri *tri;   // Triangles
    int num;    // Number of triangles
} sift3d_mesh;

/* Struct defining a 3D SIFT descriptor */
typedef struct {
    sift3d_hist hists[DESC_NUM_TOTAL_HIST]; // Array of orientation histograms
    double xd, yd, zd, sd;  // sub-pixel [x, y, z], absolute scale
} sift3d_descriptor;

/* Struct to hold SIFT3D descriptors */
typedef struct {
    sift3d_descriptor *buf;
    size_t num;
    int nx, ny, nz;         // Image dimensions
} sift3d_descriptor_store;

/* Struct to hold all parameters and internal data of the 
 * SIFT3D algorithms */
typedef struct {
    // Triange mesh
    sift3d_mesh mesh;
    // Filters for computing the GSS pyramid
    sift3d_gss_filters gss;
    // Gaussian pyramid
    sift3d_pyramid gpyr;
    // DoG pyramid
    sift3d_pyramid dog;
    // Image to process
    sift3d_image im;
    // Parameters
    double peak_thresh; // Keypoint peak threshold
    double corner_thresh; // Keypoint corner threshold
    int dense_rotate; // If true, dense descriptors are rotation-invariant
} sift3d_detector;

#ifdef __cplusplus
}
#endif

#endif

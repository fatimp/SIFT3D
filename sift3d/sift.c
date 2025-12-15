/* -----------------------------------------------------------------------------
 * sift.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains all routines needed to initialize, delete, 
 * and run the sift3d_detector detector and descriptor. It also contains routines for
 * matching sift3d_detector features and drawing the results.
 * -----------------------------------------------------------------------------
 */

#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <sift.h>
#include "imtypes_private.h"
#include "immacros.h"
#include "imutil_private.h"

/* Implementation options */
//#define SIFT3D_ORI_SOLID_ANGLE_WEIGHT // Weight bins by solid angle
//#define SIFT3D_MATCH_MAX_DIST 0.3 // Maximum distance between matching features 
//#define CUBOID_EXTREMA // Search for extrema in a cuboid region

/* Internal return codes */
#define REJECT 1

/* Default sift3d_detector parameters. These may be overriden by 
 * the calling appropriate functions. */
const double peak_thresh_default = 0.03; // DoG peak threshold
const int num_kp_levels_default = 3; // Number of levels per octave in which keypoints are found
const double corner_thresh_default = 0.4; // Minimum corner score
const double sigma_n_default = 1.15; // Nominal scale of input data
const double sigma0_default = 1.6; // Scale of the base octave

/* Internal parameters */
const double max_eig_ratio =  0.90; // Maximum ratio of eigenvalue magnitudes
const double ori_grad_thresh = 1E-10;   // Minimum norm of average gradient
const double bary_eps = FLT_EPSILON * 1E1;  // Error tolerance for barycentric coordinates
const double ori_sig_fctr = 1.5;        // Ratio of window parameter to keypoint scale
const double ori_rad_fctr =  3.0; // Ratio of window radius to parameter
const double desc_sig_fctr = 7.071067812; // See ori_sig_fctr, 5 * sqrt(2)
const double desc_rad_fctr = 2.0;  // See ori_rad_fctr
const double trunc_thresh = 0.2f * 128.0f / DESC_NUMEL; // Descriptor truncation threshold

/* Internal math constants */
const double gr = 1.6180339887; // Golden ratio

/* Get the index of bin j from triangle i */
#define MESH_GET_IDX(mesh, i, j)                \
    ((mesh)->tri[i].idx[j])

/* Get bin j from triangle i */
#define MESH_HIST_GET(mesh, hist, i, j)         \
    ((hist)->bins[MESH_GET_IDX(mesh, i, j)])

/* Clamp out of bounds polar accesses to the first or last element.
 * Note that the polar histogram is NOT circular. */
#define HIST_GET_PO(hist, a, p)                             \
    ((p) < 0 ?                                              \
     HIST_GET(hist, ((a) + NBINS_AZ / 2) % NBINS_AZ, 1) :   \
     (p) >= NBINS_PO ?                                      \
     HIST_GET(hist, ((a) + NBINS_AZ / 2) % NBINS_AZ,        \
              NBINS_PO - 1) :                               \
     HIST_GET(hist, a, p))

/* Convert out of bounds azimuthal accesses circularly, e.g. -1 goes
 * to NBINS_AZ - 1, NBINS_AZ goes to 0. This algorithm does not work
 * if the indices wrap around more than once. */
#define HIST_GET_AZ(hist, a, p)                         \
    HIST_GET_PO(hist, ((a) + NBINS_AZ) % NBINS_AZ, p)

/* Loop through a spherical image region. im and [x, y, z] are defined as
 * above. vcenter is a pointer to a sift3d_cvec specifying the center of the window.
 * rad is the radius of the window. vdisp is a pointer to a sift3d_cvec storing
 * the displacement from the window center. sqdisp is a float storing the
 * squared Euclidean distance from the window center.
 *
 * Note that the sphere is defined in real-world coordinates, i.e. those
 * with units (1, 1, 1). Thus, rad, sq_dist, and vdisp are defined in these
 * coordinates as well. However, x, y, z, and vcenter are defined in image
 * space.
 *
 * Delimit with IM_LOOP_SPHERE_END. */
#define IM_LOOP_SPHERE_START(im, x, y, z, vcenter, rad, vdisp, sq_dist) \
    {                                                                   \
    const float uxf = (float) (im)->ux;                                 \
    const float uyf = (float) (im)->uy;                                 \
    const float uzf = (float) (im)->uz;                                 \
    const int x_start = SIFT3D_MAX(floorf((vcenter)->x - (rad) / uxf), 1); \
    const int x_end   = SIFT3D_MIN(ceilf((vcenter)->x + (rad) / uxf),   \
                                   im->nx - 2);                         \
    const int y_start = SIFT3D_MAX(floorf((vcenter)->y - (rad) / uyf), 1); \
    const int y_end   = SIFT3D_MIN(ceilf((vcenter)->y + (rad) / uyf),   \
                                   im->ny - 2);                         \
    const int z_start = SIFT3D_MAX(floorf((vcenter)->z - (rad) / uzf), 1); \
    const int z_end   = SIFT3D_MIN(ceilf((vcenter)->z + (rad) / uzf),   \
                                   im->nz - 2);                         \
    SIFT3D_IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, y_start,  \
                                 y_end, z_start, z_end)                 \
        (vdisp)->x = ((float) x - (vcenter)->x) * uxf;                  \
    (vdisp)->y = ((float) y - (vcenter)->y) * uyf;                      \
    (vdisp)->z = ((float) z - (vcenter)->z) * uzf;                      \
    (sq_dist) = SIFT3D_CVEC_L2_NORM_SQ(vdisp);                          \
    if ((sq_dist) > (rad) * (rad))                                      \
        continue;                                                       \

#define IM_LOOP_SPHERE_END SIFT3D_IM_LOOP_END }

// Loop over all bins in a gradient histogram. If ICOS_HIST is defined, p
// is not referenced
#ifdef ICOS_HIST
#define HIST_LOOP_START(a, p)                           \
    for ((a) = 0; (a) < HIST_NUMEL; (a)++) { p = p; {
#else
#define HIST_LOOP_START(a, p)                   \
    for ((p) = 0; (p) < NBINS_PO; (p)++) {      \
    for ((a) = 0; (a) < NBINS_AZ; (a)++) {
#endif

// Delimit a HIST_LOOP
#define HIST_LOOP_END }}

// Get an element from a gradient histogram. If ICOS_HIST is defined, p
// is not referenced
#ifdef ICOS_HIST
#define HIST_GET_IDX(a, p) (a)
#else
#define HIST_GET_IDX(a, p) ((a) + (p) * NBINS_AZ)
#endif
#define HIST_GET(hist, a, p) ((hist)->bins[HIST_GET_IDX(a, p)])

// Get a column index in the matrix representation of a 
// sift3d_descriptor_store struct
#define DESC_MAT_GET_COL(hist_idx, a, p)                        \
    (((hist_idx) * HIST_NUMEL) + HIST_GET_IDX(a, p) + IM_NDIMS)

// As SIFT3D_IM_GET_GRAD, but with physical units (1, 1, 1)
#define IM_GET_GRAD_ISO(im, x, y, z, c, vd) {   \
        SIFT3D_IM_GET_GRAD(im, x, y, z, c, vd); \
        (vd)->x *=  1.0f / (float) (im)->ux;    \
        (vd)->y *= 1.0f / (float) (im)->uy;     \
        (vd)->z *= 1.0f / (float) (im)->uz;     \
    }

/* Initialize geometry tables. */
static int init_geometry(sift3d_detector *sift3d) {

    sift3d_mat_rm V, F;
    sift3d_cvec temp1, temp2, temp3, n;
    float mag;
    int i, j;

    sift3d_mesh * const mesh = &sift3d->mesh;

    /* Verices of a regular icosahedron inscribed in the unit sphere. */
    const float vert[] = {  0,  1,  gr,
                            0, -1,  gr,
                            0,  1, -gr,
                            0, -1, -gr,
                            1,  gr,  0,
                            -1,  gr,  0,
                            1, -gr,  0,
                            -1, -gr,  0,
                            gr,   0,  1,
                            -gr,   0,  1,
                            gr,   0, -1, 
                            -gr,   0, -1 }; 

    /* Vertex triplets forming the faces of the icosahedron. */
    const float faces[] = {0, 1, 8,
                           0, 8, 4,
                           0, 4, 5,
                           0, 5, 9,
                           0, 9, 1,
                           1, 6, 8,
                           8, 6, 10,
                           8, 10, 4,
                           4, 10, 2,
                           4, 2, 5,
                           5, 2, 11,
                           5, 11, 9,
                           9, 11, 7,
                           9, 7, 1,
                           1, 7, 6,
                           3, 6, 7,
                           3, 7, 11,
                           3, 11, 2,
                           3, 2, 10,
                           3, 10, 6};

    // Initialize matrices
    if (init_Mat_rm_p(&V, vert, ICOS_NVERT, 3, SIFT3D_FLOAT, 
                      SIFT3D_FALSE) ||
        init_Mat_rm_p(&F, faces, ICOS_NFACES, 3, SIFT3D_FLOAT, 
                      SIFT3D_FALSE))
        return SIFT3D_FAILURE;
                
    // Initialize triangle memory
    init_Mesh(mesh);
    if ((mesh->tri = (sift3d_tri *) SIFT3D_safe_realloc(mesh->tri, 
                                                 ICOS_NFACES * sizeof(sift3d_tri))) == NULL)
        return SIFT3D_FAILURE;
 
    // Populate the triangle struct for each face
    for (i = 0; i < ICOS_NFACES; i++) {

        sift3d_tri * const tri = mesh->tri + i;    
        sift3d_cvec * const v = tri->v;

        // Initialize the vertices
        for (j = 0; j < 3; j++) {

            const float mag_expected = sqrt(1 + gr * gr);
            int * const idx = tri->idx + j;

            *idx = SIFT3D_MAT_RM_GET(&F, i, j, float);

            // Initialize the vector
            v[j].x = SIFT3D_MAT_RM_GET(&V, *idx, 0, float);
            v[j].y = SIFT3D_MAT_RM_GET(&V, *idx, 1, float);
            v[j].z = SIFT3D_MAT_RM_GET(&V, *idx, 2, float);

            // Normalize to unit length
            mag = SIFT3D_CVEC_L2_NORM(v + j);
            assert(fabsf(mag - mag_expected) < 1E-10);
            SIFT3D_CVEC_SCALE(v + j, 1.0f / mag);
        }

        // Compute the normal vector at v[0] as  (V2 - V1) X (V1 - V0)
        SIFT3D_CVEC_OP(v + 2, v + 1, -, &temp1);
        SIFT3D_CVEC_OP(v + 1, v, -, &temp2);
        SIFT3D_CVEC_CROSS(&temp1, &temp2, &n);

        // Ensure this vector is facing outward from the origin
        if (SIFT3D_CVEC_DOT(&n, v) < 0) {
            // Swap two vertices
            temp1 = v[0];
            v[0] = v[1];
            v[1] = temp1;

            // Compute the normal again
            SIFT3D_CVEC_OP(v + 2, v + 1, -, &temp1);
            SIFT3D_CVEC_OP(v + 1, v, -, &temp2);
            SIFT3D_CVEC_CROSS(&temp1, &temp2, &n);
        }
        assert(SIFT3D_CVEC_DOT(&n, v) >= 0);

        // Ensure the triangle is equilateral
        SIFT3D_CVEC_OP(v + 2, v, -, &temp3);
        assert(fabsf(SIFT3D_CVEC_L2_NORM(&temp1) - 
                     SIFT3D_CVEC_L2_NORM(&temp2)) < 1E-10);
        assert(fabsf(SIFT3D_CVEC_L2_NORM(&temp1) - 
                     SIFT3D_CVEC_L2_NORM(&temp3)) < 1E-10);
    }   
    
    return SIFT3D_SUCCESS;
}

/* Convert Cartesian coordinates to barycentric. bary is set to all zeros if
 * the problem is unstable. 
 *
 * The output value k is the constant by which the ray is multiplied to
 * intersect the supporting plane of the triangle.
 *
 * This code uses the Moller-Trumbore algorithm. */
static int cart2bary(const sift3d_cvec * const cart, const sift3d_tri * const tri, 
                     sift3d_cvec * const bary, float * const k) {

    sift3d_cvec e1, e2, t, p, q;
    float det, det_inv;

    const sift3d_cvec * const v = tri->v;

    SIFT3D_CVEC_OP(v + 1, v, -, &e1);
    SIFT3D_CVEC_OP(v + 2, v, -, &e2);
    SIFT3D_CVEC_CROSS(cart, &e2, &p);
    det = SIFT3D_CVEC_DOT(&e1, &p);

    // Reject unstable points
    if (fabsf(det) < bary_eps) {
        return SIFT3D_FAILURE;
    }

    det_inv = 1.0f / det;

    t = v[0];
    SIFT3D_CVEC_SCALE(&t, -1.0f);   

    SIFT3D_CVEC_CROSS(&t, &e1, &q);

    bary->y = det_inv * SIFT3D_CVEC_DOT(&t, &p);    
    bary->z = det_inv * SIFT3D_CVEC_DOT(cart, &q);
    bary->x = 1.0f - bary->y - bary->z;

    *k = SIFT3D_CVEC_DOT(&e2, &q) * det_inv;

#ifndef NDEBUG
    sift3d_cvec temp1, temp2, temp3;
    double residual;

    if (isnan(bary->x) || isnan(bary->y) || isnan(bary->z)) {
        printf("cart2bary: invalid bary (%f, %f, %f)\n", bary->x, 
               bary->y, bary->z);
        //exit(1);
    }

    // Verify k * c = bary->x * v1 + bary->y * v2 + bary->z * v3
    temp1 = v[0];
    temp2 = v[1];
    temp3 = v[2];
    SIFT3D_CVEC_SCALE(&temp1, bary->x);
    SIFT3D_CVEC_SCALE(&temp2, bary->y); 
    SIFT3D_CVEC_SCALE(&temp3, bary->z); 
    SIFT3D_CVEC_OP(&temp1, &temp2, +, &temp1);
    SIFT3D_CVEC_OP(&temp1, &temp3, +, &temp1);
    SIFT3D_CVEC_SCALE(&temp1, 1.0f / *k);
    SIFT3D_CVEC_OP(&temp1, cart, -, &temp1);
    residual = SIFT3D_CVEC_L2_NORM(&temp1);
    if (residual > bary_eps) {
        printf("cart2bary: residual: %f\n", residual);
        exit(1);
    }
#endif
    return SIFT3D_SUCCESS;
}

/* Initialize a sift3d_keypoint_store for first use.
 * This does not need to be called to reuse the store
 * for a new image. */
static void init_Keypoint_store(sift3d_keypoint_store *const kp) {
    init_Slab(&kp->slab);
    kp->buf = (sift3d_keypoint *) kp->slab.buf;
}

/* Initialize a sift3d_keypoint struct for use. This sets up the internal pointers,
 * and nothing else. If called on a valid sift3d_keypoint struct, it has no effect. */
static int init_Keypoint(sift3d_keypoint *const key) {
    // Initialize the orientation matrix with static memory
    return init_Mat_rm_p(&key->R, key->r_data, IM_NDIMS, IM_NDIMS, 
                         SIFT3D_FLOAT, SIFT3D_FALSE);
}

/* Make room for at least num sift3d_keypoint structs in kp. 
 * 
 * Note: This function must re-initialize some internal data if it was moved. 
 * This does not affect the end user, but it affects the implementation of 
 * init_Keypoint. */
static int resize_Keypoint_store(sift3d_keypoint_store *const kp, const size_t num) {

    void *const buf_old = kp->slab.buf;

    // Resize the internal memory
    SIFT3D_RESIZE_SLAB(&kp->slab, num, sizeof(sift3d_keypoint));
    kp->buf = kp->slab.buf; 

    // If the size has changed, re-initialize the keypoints
    if (buf_old != kp->slab.buf) { 
        int i; 
        for (i = 0; i < kp->slab.num; i++) { 
            sift3d_keypoint *const key = kp->buf + i; 
            if (init_Keypoint(key)) 
                return SIFT3D_FAILURE; 
        } 
    } 

    return SIFT3D_SUCCESS;
}

/* Copy one sift3d_keypoint struct into another. */
static int copy_Keypoint(const sift3d_keypoint *const src, sift3d_keypoint *const dst) {

    // Copy the shallow data 
    dst->xd = src->xd;
    dst->yd = src->yd;
    dst->zd = src->zd;
    dst->sd = src->sd;
    dst->o = src->o;
    dst->s = src->s;

    // Copy the orienation matrix
    return copy_Mat_rm(&src->R, &dst->R);
}

/* Free all memory associated with a Keypoint_store. kp cannot be
 * used after calling this function, unless re-initialized. */
static void cleanup_Keypoint_store(sift3d_keypoint_store *const kp) {
    cleanup_Slab(&kp->slab);
}

/* Initialize a SIFT_Descriptor_store for first use.
 * This does not need to be called to reuse the store
 * for a new image. */
static void init_SIFT3D_Descriptor_store(sift3d_descriptor_store *const desc) {
    desc->buf = NULL;
}

/* Free all memory associated with a SIFT3D_Descriptor_store. desc
 * cannot be used after calling this function, unless re-initialized. */
static void cleanup_SIFT3D_Descriptor_store(sift3d_descriptor_store *const desc) {
    free(desc->buf);
}

/* Resize a sift3d_descriptor_store to hold n descriptors. Must be initialized
 * prior to calling this function. num must be positive.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
static int resize_SIFT3D_Descriptor_store(sift3d_descriptor_store *const desc,
                                          const int num) {

    if (num < 1) {
        SIFT3D_ERR("resize_SIFT3D_Descriptor_store: invalid size: %d",
                   num);
        return SIFT3D_FAILURE;
    }

    if ((desc->buf = SIFT3D_safe_realloc(desc->buf, num * sizeof(sift3d_descriptor))) == NULL)
        return SIFT3D_FAILURE;

    desc->num = num;
    return SIFT3D_SUCCESS;
}

/* Resize a sift3d_detector struct, allocating temporary storage and recompiling the 
 * filters. Does nothing unless set_im_SIFT3D was previously called. */
static int resize_SIFT3D(sift3d_detector *const sift3d, const int num_kp_levels) {

    int num_octaves; 

    const sift3d_image *const im = &sift3d->im;
    sift3d_pyramid *const gpyr = &sift3d->gpyr;
    sift3d_pyramid *const dog = &sift3d->dog;
    const unsigned int num_dog_levels = num_kp_levels + 2;
    const unsigned int num_gpyr_levels = num_dog_levels + 1;
    const int first_octave = 0;
    const int first_level = -1;

    // Compute the meximum allowed number of octaves
    if (im->data != NULL) {
        // The minimum size of a pyramid level is 8 in any dimension
        const int last_octave = 
            (int) log2((double) SIFT3D_MIN(SIFT3D_MIN(im->nx, im->ny), 
                                           im->nz)) - 3 - first_octave;

        // Verify octave parameters
        if (last_octave < first_octave) {
            SIFT3D_ERR("resize_SIFT3D: input image is too small: "
                       "must have at least 8 voxels in each "
                       "dimension \n");
            return SIFT3D_FAILURE;
        }

        num_octaves = last_octave - first_octave + 1;
    } else {
        num_octaves = 0;
    }

    // Resize the pyramid
    if (resize_Pyramid(im, first_level, num_kp_levels,
                       num_gpyr_levels, first_octave, num_octaves, gpyr) ||
        resize_Pyramid(im, first_level, num_kp_levels, 
                       num_dog_levels, first_octave, num_octaves, dog))
        return SIFT3D_FAILURE;

    // Do nothing more if we have no image
    if (im->data == NULL)
        return SIFT3D_SUCCESS;

    // Compute the Gaussian filters
    if (make_gss(&sift3d->gss, &sift3d->gpyr))
        return SIFT3D_FAILURE;

    return SIFT3D_SUCCESS;
}

/* Helper function to set the scale parameters for a sift3d_detector struct. */
static int set_scales_SIFT3D(sift3d_detector *const sift3d, const double sigma0,
                             const double sigma_n) {

    sift3d_pyramid *const gpyr = &sift3d->gpyr;
    sift3d_pyramid *const dog = &sift3d->dog;
    sift3d_gss_filters *const gss = &sift3d->gss;

    // Set the scales for the GSS and DOG pyramids
    if (set_scales_Pyramid(sigma0, sigma_n, gpyr) ||
        set_scales_Pyramid(sigma0, sigma_n, dog))
        return SIFT3D_FAILURE;

    // Do nothing more if we have no image
    if (sift3d->im.data == NULL)
        return SIFT3D_SUCCESS;

    // Recompute the filters
    return make_gss(gss, gpyr);
}

/* Sets the peak threshold, checking that it is in the interval (0, inf) */
int sift3d_detector_set_peak_thresh(sift3d_detector *const sift3d,
                                    const double peak_thresh) {
    if (peak_thresh <= 00 || peak_thresh > 1) {
        SIFT3D_ERR("sift3d_detector peak_thresh must be in the interval (0, 1]. "
                   "Provided: %f \n", peak_thresh);
        return SIFT3D_FAILURE;
    }

    sift3d->peak_thresh = peak_thresh;
    return SIFT3D_SUCCESS;
}

/* Sets the corner threshold, checking that it is in the interval [0, 1]. */
int sift3d_detector_set_corner_thresh(sift3d_detector *const sift3d,
                                      const double corner_thresh) {

    if (corner_thresh < 0.0 || corner_thresh > 1.0) {
        SIFT3D_ERR("sift3d_detector corner_thresh must be in the interval "
                   "[0, 1]. Provided: %f \n", corner_thresh);
        return SIFT3D_FAILURE;
    }

    sift3d->corner_thresh = corner_thresh;
    return SIFT3D_SUCCESS;
}

/* Sets the number of levels per octave. This function will resize the
 * internal data. */
int sift3d_detector_set_num_kp_levels(sift3d_detector *const sift3d,
                                      const unsigned int num_kp_levels) {

    const sift3d_pyramid *const gpyr = &sift3d->gpyr;

    return resize_SIFT3D(sift3d, num_kp_levels);
}

/* Sets the nominal scale parameter of the input data, checking that it is 
 * nonnegative. */
int sift3d_detector_set_sigma_n(sift3d_detector *const sift3d,
                                const double sigma_n) {

    const double sigma0 = sift3d->gpyr.sigma0;

    if (sigma_n < 0.0) {
        SIFT3D_ERR("sift3d_detector sigma_n must be nonnegative. Provided: "
                   "%f \n", sigma_n);
        return SIFT3D_FAILURE;
    }

    return set_scales_SIFT3D(sift3d, sigma0, sigma_n);
}

/* Sets the scale parameter of the first level of octave 0, checking that it
 * is nonnegative. */
int sift3d_detector_set_sigma0(sift3d_detector *const sift3d,
                               const double sigma0) {

    const double sigma_n = sift3d->gpyr.sigma_n;

    if (sigma0 < 0.0) {
        SIFT3D_ERR("sift3d_detector sigma0 must be nonnegative. Provided: "
                   "%f \n", sigma0);
        return SIFT3D_FAILURE; 
    } 

    return set_scales_SIFT3D(sift3d, sigma0, sigma_n);
}

/* Initialize a sift3d_detector struct with the default parameters. */
static int init_SIFT3D(sift3d_detector *sift3d) {

    sift3d_pyramid *const dog = &sift3d->dog;
    sift3d_pyramid *const gpyr = &sift3d->gpyr;
    sift3d_gss_filters *const gss = &sift3d->gss;

    // Initialize to defaults
    const double peak_thresh = peak_thresh_default;
    const double corner_thresh = corner_thresh_default;
    const int num_kp_levels = num_kp_levels_default;
    const double sigma_n = sigma_n_default;
    const double sigma0 = sigma0_default;
    const int dense_rotate = SIFT3D_FALSE;

    // First-time pyramid initialization
    init_Pyramid(dog);
    init_Pyramid(gpyr);

    // First-time filter initialization
    init_GSS_filters(gss);

    // Intialize the geometry tables
    if (init_geometry(sift3d))
        return SIFT3D_FAILURE;

    // Initialize the image data
    init_im(&sift3d->im);

    // Save data
    dog->first_level = gpyr->first_level = -1;
    sift3d->dense_rotate = dense_rotate;
    if (sift3d_detector_set_sigma_n(sift3d, sigma_n) ||
        sift3d_detector_set_sigma0(sift3d, sigma0) ||
        sift3d_detector_set_peak_thresh(sift3d, peak_thresh) ||
        sift3d_detector_set_corner_thresh(sift3d, corner_thresh) ||
        sift3d_detector_set_num_kp_levels(sift3d, num_kp_levels))
        return SIFT3D_FAILURE;

    return SIFT3D_SUCCESS;
}

/* Free all memory associated with a sift3d_detector struct. sift3d cannot be reused
 * unless it is reinitialized. */
static void cleanup_SIFT3D(sift3d_detector *const sift3d) {

    // Clean up the image copy
    im_free(&sift3d->im);

    // Clean up the pyramids
    cleanup_Pyramid(&sift3d->gpyr);
    cleanup_Pyramid(&sift3d->dog);

    // Clean up the GSS filters
    cleanup_GSS_filters(&sift3d->gss);

    // Clean up the triangle mesh 
    cleanup_Mesh(&sift3d->mesh);
}

/* Helper routine to begin processing a new image. If the dimensions differ
 * from the last one, this function resizes the sift3d_detector struct. */
static int set_im_SIFT3D(sift3d_detector *const sift3d, const sift3d_image *const im) {

    int dims_old[IM_NDIMS];
    int i;

    const float *const data_old = sift3d->im.data;
    const sift3d_pyramid *const gpyr = &sift3d->gpyr;
    const int first_octave = sift3d->gpyr.first_octave;
    const int num_kp_levels = gpyr->num_kp_levels;

    // Make a temporary copy the previous image dimensions
    for (i = 0; i < IM_NDIMS; i++) {
        dims_old[i] = SIFT3D_IM_GET_DIMS(&sift3d->im)[i];
    }

    // Make a copy of the input image
    if (im_copy_data(im, &sift3d->im))
        return SIFT3D_FAILURE;

    // Scale the input image to [-1, 1]
    im_scale(&sift3d->im);

    // Resize the internal data, if necessary
    if ((data_old == NULL || 
         memcmp(dims_old, SIFT3D_IM_GET_DIMS(&sift3d->im), 
                IM_NDIMS * sizeof(int))) &&
        resize_SIFT3D(sift3d, num_kp_levels))
        return SIFT3D_FAILURE;

    return SIFT3D_SUCCESS;
}

/* Build the GSS pyramid on a single CPU thread */
static int build_gpyr(sift3d_detector *sift3d) {

    const sift3d_image *prev;
    sift3d_sep_fir_filter *f;
    sift3d_image *cur;
    int o, s;

    sift3d_pyramid *const gpyr = &sift3d->gpyr;
    const sift3d_gss_filters *const gss = &sift3d->gss;
    const int s_start = gpyr->first_level + 1;
    const int s_end = SIFT3D_PYR_LAST_LEVEL(gpyr);
    const int o_start = gpyr->first_octave;
    const int o_end = SIFT3D_PYR_LAST_OCTAVE(gpyr);
    const double unit = 1.0;

    // Build the first image
    cur = SIFT3D_PYR_IM_GET(gpyr, o_start, s_start - 1);
    prev = &sift3d->im;

    f = (sift3d_sep_fir_filter *) &gss->first_gauss.f;
    if (apply_Sep_FIR_filter(prev, cur, f, unit))
        return SIFT3D_FAILURE;

    // Build the rest of the pyramid
    SIFT3D_PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end)
        cur = SIFT3D_PYR_IM_GET(gpyr, o, s);
    prev = SIFT3D_PYR_IM_GET(gpyr, o, s - 1);
    f = &gss->gauss_octave[s].f;
    if (apply_Sep_FIR_filter(prev, cur, f, unit))
        return SIFT3D_FAILURE;
    SIFT3D_PYR_LOOP_SCALE_END
        // Downsample
        if (o != o_end) {

            const int downsample_level = 
                SIFT3D_MAX(s_end - 2, gpyr->first_level);

            prev = SIFT3D_PYR_IM_GET(gpyr, o, downsample_level);
            cur = SIFT3D_PYR_IM_GET(gpyr, o + 1, s_start - 1);

            assert(fabs(prev->s - cur->s) < FLT_EPSILON);

            if (im_downsample_2x(prev, cur))
                return SIFT3D_FAILURE;

        }
    SIFT3D_PYR_LOOP_OCTAVE_END

    return SIFT3D_SUCCESS;
}

static int build_dog(sift3d_detector *sift3d) {

    sift3d_image *gpyr_cur, *gpyr_next, *dog_level;
    int o, s;

    sift3d_pyramid *const dog = &sift3d->dog;
    sift3d_pyramid *const gpyr = &sift3d->gpyr;

    SIFT3D_PYR_LOOP_START(dog, o, s)
        gpyr_cur = SIFT3D_PYR_IM_GET(gpyr, o, s);
    gpyr_next = SIFT3D_PYR_IM_GET(gpyr, o, s + 1);          
    dog_level = SIFT3D_PYR_IM_GET(dog, o, s);
        
    if (im_subtract(gpyr_cur, gpyr_next, 
                    dog_level))
        return SIFT3D_FAILURE;
    SIFT3D_PYR_LOOP_END

        return SIFT3D_SUCCESS;
}

/* Detect local extrema */
static int detect_extrema(sift3d_detector *sift3d, sift3d_keypoint_store *kp) {

    sift3d_image *cur, *prev, *next;
    sift3d_keypoint *key;
    float pcur, peak_thresh;
    int o, s, x, y, z, x_start, x_end, y_start, y_end, z_start, z_end, num;

    const sift3d_pyramid *const dog = &sift3d->dog;
    const int o_start = dog->first_octave;
    const int o_end = SIFT3D_PYR_LAST_OCTAVE(dog);
    const int s_start = dog->first_level + 1;
    const int s_end = SIFT3D_PYR_LAST_LEVEL(dog) - 1;

    // Verify the inputs
    if (dog->num_levels < 3) {
        printf("detect_extrema: Requires at least 3 levels per octave, "
               "provided only %d \n", dog->num_levels);
        return SIFT3D_FAILURE;
    }

    // Initialize dimensions of keypoint store
    cur = SIFT3D_PYR_IM_GET(dog, o_start, s_start);
    kp->nx = cur->nx;
    kp->ny = cur->ny;
    kp->nz = cur->nz;

#define CMP_CUBE(im, x, y, z, CMP, IGNORESELF, val) (                   \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) - 1, (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) - 1, (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) + 1, (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) + 1, (z) - 1, 0) && \
        ((val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),    (z), 0   ) || \
         IGNORESELF) &&                                                 \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) - 1, (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) - 1, (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) + 1, (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) + 1, (z), 0    ) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) - 1, (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) - 1, (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) + 1, (z) + 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) + 1, (z) + 1, 0) )
#ifdef CUBOID_EXTREMA
#define CMP_PREV(im, x, y, z, CMP, val)             \
    CMP_CUBE(im, x, y, z, CMP, SIFT3D_FALSE, val)
#define CMP_CUR(im, x, y, z, CMP, val)              \
    CMP_CUBE(im, x, y, z, CMP, SIFT3D_TRUE, val)
#define CMP_NEXT(im, x, y, z, CMP, val)             \
    CMP_CUBE(im, x, y, z, CMP, SIFT3D_FALSE, val)
#else
#define CMP_PREV(im, x, y, z, CMP, val) (                       \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x), (y), (z), 0)    \
        )
#define CMP_CUR(im, x, y, z, CMP, val) (                                \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z), 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z), 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z), 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z), 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) - 1, 0) && \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) + 1, 0) \
        )
#define CMP_NEXT(im, x, y, z, CMP, val)         \
    CMP_PREV(im, x, y, z, CMP, val)
#endif

    num = 0;
    SIFT3D_PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end)  

        // Select current and neighboring levels
        prev = SIFT3D_PYR_IM_GET(dog, o, s - 1);
    cur = SIFT3D_PYR_IM_GET(dog, o, s);
    next = SIFT3D_PYR_IM_GET(dog, o, s + 1);
    peak_thresh = sift3d->peak_thresh;

    // Loop through all non-boundary pixels
    x_start = y_start = z_start = 1;
    x_end = cur->nx - 2;
    y_end = cur->ny - 2;
    z_end = cur->nz - 2;
    SIFT3D_IM_LOOP_LIMITED_START(cur, x, y, z, x_start, x_end, y_start,
                                 y_end, z_start, z_end)
        // Sample the center value
        pcur = SIFT3D_IM_GET_VOX(cur, x, y, z, 0);

    // Apply the peak threshold
    if ((pcur > peak_thresh || pcur < -peak_thresh) &&
        // Compare to the neighbors
        ((CMP_PREV(prev, x, y, z, >, pcur) &&
          CMP_CUR(cur, x, y, z, >, pcur) &&
          CMP_NEXT(next, x, y, z, >, pcur)) ||
         (CMP_PREV(prev, x, y, z, <, pcur) &&
          CMP_CUR(cur, x, y, z, <, pcur) &&
          CMP_NEXT(next, x, y, z, <, pcur)))) {

        // Add a keypoint candidate
        num++;
        if (resize_Keypoint_store(kp, num))
            return SIFT3D_FAILURE;
        key = kp->buf + num - 1;
        if (init_Keypoint(key))
            return SIFT3D_FAILURE;
        key->o = o;
        key->s = s;
        key->sd = cur->s;
        key->xd = (double) x;
        key->yd = (double) y;
        key->zd = (double) z;
        key->strength = fabsf (pcur);
    }
    SIFT3D_IM_LOOP_END
        SIFT3D_PYR_LOOP_END
#undef CMP_NEIGHBORS

        return SIFT3D_SUCCESS;
}

/* Bin a Cartesian gradient into Spherical gradient bins */
SIFT3D_IGNORE_UNUSED
static int Cvec_to_sbins(const sift3d_cvec * const vd, sift3d_svec * const bins) {

    // Convert to spherical coordinates
    SIFT3D_CVEC_TO_SVEC(vd, bins);
    //FIXME: Is this needed? SIFT3D_CVEC_TO_SVEC cannot divide by zero
    if (bins->mag < FLT_EPSILON * 1E2)
        return SIFT3D_FAILURE;

    // Compute bins
    bins->az *= (float) NBINS_AZ / SIFT3D_AZ_MAX_F; 
    bins->po *= (float) NBINS_PO / SIFT3D_PO_MAX_F;

    assert(bins->az < NBINS_AZ);
    assert(bins->po < NBINS_PO);

    return SIFT3D_SUCCESS;
}

/* Refine a gradient histogram with optional operations,
 * such as solid angle weighting. */
static void refine_Hist(sift3d_hist *hist) {

#ifndef ICOS_HIST

#ifdef SIFT3D_ORI_SOLID_ANGLE_WEIGHT
    {
        float po;
        int a, p;
        // TODO: could accelerate this with a lookup table      

        // Weight by the solid angle of the bins, ignoring constants
        HIST_LOOP_START(a, p)
            po = p * po_max_f / NBINS_PO;
        HIST_GET(hist, a, p) /= cosf(po) - cosf(po + po_max_f / NBINS_PO);
        HIST_LOOP_END
            }
#endif

#endif

}

/* Assign an orientation to a point in an image.
 *
 * Parameters:
 *   -im: The image data.
 *   -vcenter: The center of the window, in image space.
 *   -sigma: The scale parameter. The width of the window is a constant
 *      multiple of this.
 *   -R: The place to write the rotation matrix.
 */
static int assign_eig_ori(const sift3d_image *const im, const sift3d_cvec *const vcenter,
                          const double sigma, sift3d_mat_rm *const R, 
                          double *const conf) {
    sift3d_cvec v[2];
    sift3d_mat_rm A, L, Q;
    sift3d_cvec vd_win, vdisp, vr;
    double d, cos_ang, abs_cos_ang, corner_score;
    float weight, sq_dist, sgn;
    int i, x, y, z, m;
  
    const double win_radius = sigma * ori_rad_fctr; 

    // Verify inputs
    if (!SIFT3D_IM_CONTAINS_CVEC(im, vcenter)) {
        SIFT3D_ERR("assign_eig_ori: vcenter (%f, %f, %f) lies "
                   "outside the boundaries of im [%d x %d x %d] \n", 
                   vcenter->x, vcenter->y, vcenter->z, im->nx, im->ny, im->nz);
        return SIFT3D_FAILURE;
    }
    if (sigma < 0) {
        SIFT3D_ERR("assign_eig_ori: invalid sigma: %f \n", sigma);
        return SIFT3D_FAILURE;
    }

    // Initialize the intermediates
    if (init_Mat_rm(&A, 3, 3, SIFT3D_DOUBLE, SIFT3D_TRUE))
        return SIFT3D_FAILURE;
    if (init_Mat_rm(&L, 0, 0, SIFT3D_DOUBLE, SIFT3D_TRUE) ||
        init_Mat_rm(&Q, 0, 0, SIFT3D_DOUBLE, SIFT3D_TRUE))
        goto eig_ori_fail;

    // Resize the output
    R->num_rows = R->num_cols = IM_NDIMS;
    R->type = SIFT3D_FLOAT;
    if (resize_Mat_rm(R))
        goto eig_ori_fail;

    // Form the structure tensor and window gradient
    vd_win.x = 0.0f;
    vd_win.y = 0.0f;
    vd_win.z = 0.0f;
    IM_LOOP_SPHERE_START(im, x, y, z, vcenter, win_radius, &vdisp, sq_dist)

        sift3d_cvec vd;

    // Compute Gaussian weighting, ignoring the constant factor
    weight = expf(-0.5 * sq_dist / (sigma * sigma));        

    // Get the gradient 
    IM_GET_GRAD_ISO(im, x, y, z, 0, &vd);

    // Update the structure tensor
    SIFT3D_MAT_RM_GET(&A, 0, 0, double) += (double) vd.x * vd.x * weight;
    SIFT3D_MAT_RM_GET(&A, 0, 1, double) += (double) vd.x * vd.y * weight;
    SIFT3D_MAT_RM_GET(&A, 0, 2, double) += (double) vd.x * vd.z * weight;
    SIFT3D_MAT_RM_GET(&A, 1, 1, double) += (double) vd.y * vd.y * weight;
    SIFT3D_MAT_RM_GET(&A, 1, 2, double) += (double) vd.y * vd.z * weight;
    SIFT3D_MAT_RM_GET(&A, 2, 2, double) += (double) vd.z * vd.z * weight;

    // Update the window gradient
    SIFT3D_CVEC_SCALE(&vd, weight);
    SIFT3D_CVEC_OP(&vd_win, &vd, +, &vd_win);

    IM_LOOP_SPHERE_END

        // Fill in the remaining elements
        SIFT3D_MAT_RM_GET(&A, 1, 0, double) = SIFT3D_MAT_RM_GET(&A, 0, 1, double);
    SIFT3D_MAT_RM_GET(&A, 2, 0, double) = SIFT3D_MAT_RM_GET(&A, 0, 2, double);
    SIFT3D_MAT_RM_GET(&A, 2, 1, double) = SIFT3D_MAT_RM_GET(&A, 1, 2, double);

    // Reject keypoints with weak gradient 
    if (SIFT3D_CVEC_L2_NORM_SQ(&vd_win) < (float) ori_grad_thresh) {
        goto eig_ori_reject;
    } 

    // Get the eigendecomposition
    if (eigen_Mat_rm(&A, &Q, &L))
        goto eig_ori_fail;

    // Ensure we have distinct eigenvalues
    m = L.num_rows;
    if (m != 3)
        goto eig_ori_reject;

    // Test the eigenvectors for stability
    for (i = 0; i < m - 1; i++) {
        if (fabs(SIFT3D_MAT_RM_GET(&L, i, 0, double) /
                 SIFT3D_MAT_RM_GET(&L, i + 1, 0, double)) > max_eig_ratio)
            goto eig_ori_reject;
    }

    // Assign signs to the first n - 1 vectors
    corner_score = DBL_MAX;
    for (i = 0; i < m - 1; i++) {

        const int eig_idx = m - i - 1;

        // Get an eigenvector, in descending order
        vr.x = (float) SIFT3D_MAT_RM_GET(&Q, 0, eig_idx, double);
        vr.y = (float) SIFT3D_MAT_RM_GET(&Q, 1, eig_idx, double);
        vr.z = (float) SIFT3D_MAT_RM_GET(&Q, 2, eig_idx, double);

        // Get the directional derivative
        d = SIFT3D_CVEC_DOT(&vd_win, &vr);

        // Get the cosine of the angle between the eigenvector and the gradient
        cos_ang = d / (SIFT3D_CVEC_L2_NORM(&vr) * SIFT3D_CVEC_L2_NORM(&vd_win));
        abs_cos_ang = fabs(cos_ang);

        // Compute the corner confidence score
        corner_score = SIFT3D_MIN(corner_score, abs_cos_ang);

        // Get the sign of the derivative
        sgn = d > 0.0 ? 1.0f : -1.0f;

        // Enforce positive directional derivative
        SIFT3D_CVEC_SCALE(&vr, sgn);

        // Add the vector to the rotation matrix
        SIFT3D_MAT_RM_GET(R, 0, i, float) = vr.x;
        SIFT3D_MAT_RM_GET(R, 1, i, float) = vr.y;
        SIFT3D_MAT_RM_GET(R, 2, i, float) = vr.z;

        // Save this vector for later use
        v[i] = vr;
    }

    // Take the cross product of the first two vectors
    SIFT3D_CVEC_CROSS(v, v + 1, &vr);

    // Add the last vector
    SIFT3D_MAT_RM_GET(R, 0, 2, float) = (float) vr.x;
    SIFT3D_MAT_RM_GET(R, 1, 2, float) = (float) vr.y;
    SIFT3D_MAT_RM_GET(R, 2, 2, float) = (float) vr.z;

    // Optionally write back the corner score
    if (conf != NULL)
        *conf = corner_score;

    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return SIFT3D_SUCCESS; 

eig_ori_reject:
    if (conf != NULL)
        *conf = 0.0;
    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return REJECT;

eig_ori_fail:
    if (conf != NULL)
        *conf = 0.0;
    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return SIFT3D_FAILURE;
}

/* Helper function to call assign_eig_ori, and reject keypoints with
 * confidence below the parameter "thresh." All other parameters are the same.
 * All return values are the same, except REJECT is returned if 
 * conf < thresh. */
static int assign_orientation_thresh(const sift3d_image *const im, 
                                     const sift3d_cvec *const vcenter, const double sigma, const double thresh,
                                     sift3d_mat_rm *const R) {

    double conf;
    int ret;

    ret = assign_eig_ori(im, vcenter, sigma, R, &conf);

    return ret == SIFT3D_SUCCESS ? 
        (conf < thresh ? REJECT : SIFT3D_SUCCESS) : ret;
}

/* Assign rotation matrices to the keypoints. 
 * 
 * Note that this stage will modify kp, likely
 * rejecting some keypoints as orientationally
 * unstable. */
static int assign_orientations(sift3d_detector *const sift3d,
                               sift3d_keypoint_store *const kp) {
    sift3d_keypoint *kp_pos;
    size_t num;
    int i, err;

    // Iterate over the keypoints
    err = SIFT3D_SUCCESS;
#pragma omp parallel for
    for (i = 0; i < kp->slab.num; i++) {

        sift3d_keypoint *const key = kp->buf + i;
        const sift3d_image *const level =
            SIFT3D_PYR_IM_GET(&sift3d->gpyr, key->o, key->s);
        sift3d_mat_rm *const R = &key->R;
        const sift3d_cvec vcenter = {key->xd, key->yd, key->zd};
        const double sigma = ori_sig_fctr * key->sd;

        // Compute dominant orientations
        assert(R->u.data_float == key->r_data);
        switch (assign_orientation_thresh(level, &vcenter, sigma,
                                          sift3d->corner_thresh, R)) {
        case SIFT3D_SUCCESS:
            // Continue processing this keypoint
            break;
        case REJECT:
            // Mark this keypoint as invalid
            key->xd = key->yd = key->zd = -1.0;
            continue;
        default:
            // Any other return value is an error
            err = SIFT3D_FAILURE;
            continue;
        }
    }

    // Check for errors
    if (err) return err;

    // Rebuild the keypoint buffer in place
    kp_pos = kp->buf;
    for (i = 0; i < kp->slab.num; i++) {
        sift3d_keypoint *const key = kp->buf + i;

        // Check if the keypoint is valid
        if (key->xd < 0.0)
            continue;

        // Copy this keypoint to the next available spot
        if (copy_Keypoint(key, kp_pos))
            return SIFT3D_FAILURE;

        kp_pos++;
    }

    // Release unneeded keypoint memory
    num = kp_pos - kp->buf;
    return resize_Keypoint_store(kp, num);
}

/* Verify that keypoints kp are valid in image im. Returns SIFT3D_SUCCESS if
 * valid, SIFT3D_FAILURE otherwise. */
static int verify_keys(const sift3d_keypoint_store *const kp, const sift3d_image *const im) {

    int i;

    const int num = kp->slab.num;

    // Check the number of keypoints
    if (num < 1) {
        SIFT3D_ERR("verify_keys: invalid number of keypoints: "
                   "%d \n", num);
        return SIFT3D_FAILURE;
    }

    // Check each keypoint
    for (i = 0; i < num; i++) {

        const sift3d_keypoint *key = kp->buf + i;

        const double octave_factor = ldexp(1.0, key->o);

        if (key->xd < 0 ||
            key->yd < 0 ||
            key->zd < 0 ||
            key->xd * octave_factor >= (double) im->nx || 
            key->yd * octave_factor >= (double) im->ny || 
            key->zd * octave_factor >= (double) im->nz) {
            SIFT3D_ERR("verify_keys: keypoint %d (%f, %f, %f) "
                       "octave %d exceeds image dimensions "
                       "(%d, %d, %d) \n", i, key->xd, key->yd, key->zd,
                       key->o, im->nx, im->ny, im->nz);
            return SIFT3D_FAILURE; 
        }

        if (key->sd <= 0) {
            SIFT3D_ERR("verify_keys: keypoint %d has invalid "
                       "scale %f \n", i, key->sd);
            return SIFT3D_FAILURE;
        }
    }

    return SIFT3D_SUCCESS;
}

/* Detect keypoint locations and orientations. You must initialize
 * the sift3d_detector struct, image, and keypoint store with the appropriate
 * functions prior to calling this function. */
int sift3d_detect_keypoints(sift3d_detector *const sift3d,
                            const sift3d_image *const im,
                            sift3d_keypoint_store *const kp) {
    // Verify inputs
    if (im->nc != 1) {
        SIFT3D_ERR("SIFT3D_detect_keypoints: invalid number "
                   "of image channels: %d -- only single-channel images "
                   "are supported \n", im->nc);
        return SIFT3D_FAILURE;
    }

    // Set the image       
    if (set_im_SIFT3D(sift3d, im))
        return SIFT3D_FAILURE;

    // Build the GSS pyramid
    if (build_gpyr(sift3d))
        return SIFT3D_FAILURE;

    // Build the DoG pyramid
    if (build_dog(sift3d))
        return SIFT3D_FAILURE;

    // Detect extrema
    if (detect_extrema(sift3d, kp))
        return SIFT3D_FAILURE;

    // Assign orientations
    if (assign_orientations(sift3d, kp))
        return SIFT3D_FAILURE;

    return SIFT3D_SUCCESS;
}

/* Get the bin and barycentric coordinates of a vector in the icosahedral 
 * histogram. */
SIFT3D_IGNORE_UNUSED
static int icos_hist_bin(const sift3d_detector * const sift3d,
                         const sift3d_cvec * const x, sift3d_cvec * const bary,
                         int * const bin) { 

    float k;
    int i;

    const sift3d_mesh * const mesh = &sift3d->mesh;

    // Check for very small vectors
    if (SIFT3D_CVEC_L2_NORM_SQ(x) < bary_eps)
        return SIFT3D_FAILURE;

    // Iterate through the faces
    for (i = 0; i < ICOS_NFACES; i++) {

        const sift3d_tri * const tri = mesh->tri + i;

        // Convert to barycentric coordinates
        if (cart2bary(x, tri, bary, &k))
            continue;

        // Test for intersection
        if (bary->x < -bary_eps || bary->y < -bary_eps ||
            bary->z < -bary_eps || k < 0)
            continue;

        // Save the bin
        *bin = i;

        // No other triangles will be intersected
        return SIFT3D_SUCCESS;
    }   

    // Unreachable code
    assert(SIFT3D_FALSE);
    return SIFT3D_FAILURE;
}

/* Helper routine to interpolate over the histograms of a
 * sift3d_detector descriptor. */
static void SIFT3D_desc_acc_interp(const sift3d_detector * const sift3d, 
                                   const sift3d_cvec * const vbins, 
                                   const sift3d_cvec * const grad,
                                   sift3d_descriptor * const desc) {

    sift3d_cvec dvbins;
    sift3d_hist *hist;
    float weight;
    int dx, dy, dz, x, y, z;

#ifdef ICOS_HIST
    sift3d_cvec bary;
    float mag;
    int bin;    
#else
    sift3d_svec sbins, dsbins;
    int da, dp, a, p;
#endif

    const int y_stride = NHIST_PER_DIM;
    const int z_stride = NHIST_PER_DIM * NHIST_PER_DIM; 

    // Compute difference from integer bin values
    dvbins.x = vbins->x - floorf(vbins->x);
    dvbins.y = vbins->y - floorf(vbins->y);
    dvbins.z = vbins->z - floorf(vbins->z);

    // Compute the histogram bin
#ifdef ICOS_HIST
    const sift3d_mesh *const mesh = &sift3d->mesh;

    // Get the index of the intersecting face 
    if (icos_hist_bin(sift3d, grad, &bary, &bin))
        return;
    
    // Get the magnitude of the vector
    mag = SIFT3D_CVEC_L2_NORM(grad);

#else
    if (Cvec_to_sbins(grad, &sbins))
        return;
    dsbins.az = sbins.az - floorf(sbins.az);
    dsbins.po = sbins.po - floorf(sbins.po);
#endif
    
    for (dx = 0; dx < 2; dx++) {
        for (dy = 0; dy < 2; dy++) {
            for (dz = 0; dz < 2; dz++) {

                x = (int) vbins->x + dx;
                y = (int) vbins->y + dy;
                z = (int) vbins->z + dz;

                // Check boundaries
                if (x < 0 || x >= NHIST_PER_DIM ||
                    y < 0 || y >= NHIST_PER_DIM ||
                    z < 0 || z >= NHIST_PER_DIM)
                    continue;

                // Get the histogram
                hist = desc->hists + x + y * y_stride + 
                    z * z_stride;    

                assert(x + y * y_stride + z * z_stride < DESC_NUM_TOTAL_HIST);

                // Get the spatial interpolation weight
                weight = ((dx == 0) ? (1.0f - dvbins.x) : dvbins.x) *
                    ((dy == 0) ? (1.0f - dvbins.y) : dvbins.y) *
                    ((dz == 0) ? (1.0f - dvbins.z) : dvbins.z);

                /* Add the value into the histogram */
#ifdef ICOS_HIST
                assert(HIST_NUMEL == ICOS_NVERT);
                assert(bin >= 0 && bin < ICOS_NFACES);

                // Interpolate over three vertices
                MESH_HIST_GET(mesh, hist, bin, 0) += mag * weight * bary.x;
                MESH_HIST_GET(mesh, hist, bin, 1) += mag * weight * bary.y;
                MESH_HIST_GET(mesh, hist, bin, 2) += mag * weight * bary.z; 
#else
                // Iterate over all angles
                for (dp = 0; dp < 2; dp ++) {
                    for (da = 0; da < 2; da ++) {

                        a = ((int) sbins.az + da) % NBINS_AZ;
                        p = (int) sbins.po + dp;
                        if (p >= NBINS_PO) {
                            // See HIST_GET_PO
                            a = (a + NBINS_AZ / 2) % NBINS_AZ;
                            p = NBINS_PO - 1;
                        }
        
                        assert(a >= 0);
                        assert(a < NBINS_AZ);
                        assert(p >= 0);
                        assert(p < NBINS_PO);

                        HIST_GET(hist, a, p) += sbins.mag * weight *
                            ((da == 0) ? (1.0f - dsbins.az) : dsbins.az) *
                            ((dp == 0) ? (1.0f - dsbins.po) : dsbins.po);
                    }}
#endif
            }}}

}

/* Normalize a descriptor */
static void normalize_desc(sift3d_descriptor * const desc) {

    double norm; 
    int i, a, p;

    norm = 0.0;
    for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) { 

        const sift3d_hist *const hist = desc->hists + i;

        HIST_LOOP_START(a, p) 
            const float el = HIST_GET(hist, a, p);
        norm += (double) el * el;
        HIST_LOOP_END 
            }

    norm = sqrt(norm) + DBL_EPSILON; 

    for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {

        sift3d_hist *const hist = desc->hists + i;
        const float norm_inv = 1.0f / norm; 

        HIST_LOOP_START(a, p) 
            HIST_GET(hist, a, p) *= norm_inv; 
        HIST_LOOP_END 
            }
}

/* Set a histogram to zero */
static void hist_zero(sift3d_hist *hist) {

    int a, p;

    HIST_LOOP_START(a, p)
        HIST_GET(hist, a, p) = 0.0f;
    HIST_LOOP_END
        }

/* Helper routine to extract a single sift3d_detector descriptor */
static int extract_descrip(sift3d_detector *const sift3d, const sift3d_image *const im,
                           const sift3d_keypoint *const key, sift3d_descriptor *const desc) {

    float buf[IM_NDIMS * IM_NDIMS];
    sift3d_mat_rm Rt;
    sift3d_cvec vcenter, vim, vkp, vbins, grad, grad_rot;
    sift3d_hist *hist;
    float weight, sq_dist;
    int i, x, y, z, a, p;

    // Compute basic parameters 
    const float sigma = key->sd * desc_sig_fctr;
    const float win_radius = desc_rad_fctr * sigma;
    const float desc_half_width = win_radius / sqrt(2);
    const float desc_width = 2.0f * desc_half_width;
    const float desc_hist_width = desc_width / NHIST_PER_DIM;
    const float desc_bin_fctr = 1.0f / desc_hist_width;
    const double coord_factor = ldexp(1.0, key->o);

    // Invert the rotation matrix
    if (init_Mat_rm_p(&Rt, buf, IM_NDIMS, IM_NDIMS, SIFT3D_FLOAT, 
                      SIFT3D_FALSE) ||
        transpose_Mat_rm(&key->R, &Rt))
        return SIFT3D_FAILURE;

    // Zero the descriptor
    for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {
        hist = desc->hists + i;
        hist_zero(hist);
    }

    // Iterate over a sphere window in real-world coordinates 
    vcenter.x = key->xd;
    vcenter.y = key->yd;
    vcenter.z = key->zd;
    IM_LOOP_SPHERE_START(im, x, y, z, &vcenter, win_radius, &vim, sq_dist)

        // Rotate to keypoint space
        SIFT3D_MUL_MAT_RM_CVEC(&Rt, &vim, &vkp);        

    // Compute spatial bins
    vbins.x = (vkp.x + desc_half_width) * desc_bin_fctr;
    vbins.y = (vkp.y + desc_half_width) * desc_bin_fctr;
    vbins.z = (vkp.z + desc_half_width) * desc_bin_fctr;

    // Reject points outside the rectangular descriptor 
    if (vbins.x < 0 || vbins.y < 0 || vbins.z < 0 ||
        vbins.x >= (float) NHIST_PER_DIM ||
        vbins.y >= (float) NHIST_PER_DIM ||
        vbins.z >= (float) NHIST_PER_DIM)
        continue;

    // Take the gradient
    IM_GET_GRAD_ISO(im, x, y, z, 0, &grad);

    // Apply a Gaussian window
    weight = expf(-0.5f * sq_dist / (sigma * sigma));
    SIFT3D_CVEC_SCALE(&grad, weight);

    // Rotate the gradient to keypoint space
    SIFT3D_MUL_MAT_RM_CVEC(&Rt, &grad, &grad_rot);

    // Finally, accumulate bins by 5x linear interpolation
    SIFT3D_desc_acc_interp(sift3d, &vbins, &grad_rot, desc);
    IM_LOOP_SPHERE_END

        // Histogram refinement steps
        for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {
            refine_Hist(&desc->hists[i]);
        }

    // Normalize the descriptor
    normalize_desc(desc);

    // Truncate
    for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {
        hist = desc->hists + i;
        HIST_LOOP_START(a, p)
            HIST_GET(hist, a, p) = SIFT3D_MIN(HIST_GET(hist, a, p), 
                                              (float) trunc_thresh);
        HIST_LOOP_END
            }

    // Normalize again
    normalize_desc(desc);

    // Save the descriptor location in the original image
    // coordinates
    desc->xd = key->xd * coord_factor;
    desc->yd = key->yd * coord_factor;
    desc->zd = key->zd * coord_factor;
    desc->sd = key->sd;

    return SIFT3D_SUCCESS;
}

/* Check if the Gaussian scale-space pyramid in a sift3d_detector struct is valid. This
 * shall return SIFT3D_TRUE if the struct was initialized, and 
 * SIFT3D_detect_keypoints has been successfully called on it since 
 * initialization. 
 *
 * Note: sift3d must be initialized before calling this function. */
int sift3d_detector_has_gpyr(const sift3d_detector *const sift3d) {
    const sift3d_pyramid *const gpyr = &sift3d->gpyr;

    return gpyr->levels != NULL && gpyr->num_levels != 0 && 
        gpyr->num_octaves != 0;
}

/* Helper funciton to extract sift3d_detector descriptors from a list of keypoints and 
 * an image. Called by SIFT3D_extract_descriptors and 
 * SIFT3D_extract_raw_descriptors.
 *
 * parameters:
 *  sift3d - (initialized) struct defining the algorithm parameters
 *  gpyr - A Gaussian Scale-Space pyramid containing the image data
 *  kp - keypoint list populated by a feature detector 
 *  desc - (initialized) struct to hold the descriptors
 *  use_gpyr - see im for details */
static int do_extract_descriptors(sift3d_detector *const sift3d,
                                  const sift3d_pyramid *const gpyr,
                                  const sift3d_keypoint_store *const kp,
                                  sift3d_descriptor_store *const desc) {
    int i, ret;

    const sift3d_image *const first_level = 
        SIFT3D_PYR_IM_GET(gpyr, gpyr->first_octave, gpyr->first_level);

    const int num = kp->slab.num;

    // Initialize the metadata 
    desc->nx = first_level->nx; 
    desc->ny = first_level->ny; 
    desc->nz = first_level->nz; 

    // Resize the descriptor store
    if (resize_SIFT3D_Descriptor_store(desc, num))
        return SIFT3D_FAILURE;

    // Extract the descriptors
    ret = SIFT3D_SUCCESS;
#pragma omp parallel for
    for (i = 0; i < desc->num; i++) {

        const sift3d_keypoint *const key = kp->buf + i;
        sift3d_descriptor *const descrip = desc->buf + i;
        const sift3d_image *const level = 
            SIFT3D_PYR_IM_GET(gpyr, key->o, key->s);

        if (extract_descrip(sift3d, level, key, descrip)) {
            ret = SIFT3D_FAILURE;
        }
    }   

    return ret;
}

/* Extract sift3d_detector descriptors from a list of keypoints. Uses the Gaussian
 * scale-space pyramid from the previous call to SIFT3D_detect_keypoints on
 * this sift3d_detector struct. To extract from an image, see 
 * SIFT3D_extract_raw_descriptors. 
 *
 * Note: To check if SIFT3D_detect_keypoints has been called on this struct,
 * use SIFT3D_have_gpyr.
 *
 * Parameters:
 *  sift3d - (initialized) struct defining the algorithm parameters. Must have
 *      been used in some previous call to SIFT3D_detect_keypoints.
 *  kp - keypoint list populated by a feature detector 
 *  desc - (initialized) struct to hold the descriptors
 *
 * Return value:
 *  Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int sift3d_extract_descriptors(sift3d_detector *const sift3d, 
                               const sift3d_keypoint_store *const kp, 
                               sift3d_descriptor_store *const desc) {
    // Verify inputs
    if (verify_keys(kp, &sift3d->im))
        return SIFT3D_FAILURE;

    // Check if a Gaussian scale-space pyramid is available for processing
    if (!sift3d_detector_has_gpyr(sift3d)) {
        SIFT3D_ERR("SIFT3D_extract_descriptors: no Gaussian pyramid is "
                   "available. Make sure SIFT3D_detect_keypoints was "
                   "called prior to calling this function. \n");
        return SIFT3D_FAILURE;
    }

    // Extract features
    if (do_extract_descriptors(sift3d, &sift3d->gpyr, kp, desc))
        return SIFT3D_FAILURE;

    return SIFT3D_SUCCESS;
}

/* Convert a keypoint store to a matrix. 
 * Output format:
 *  [x1 y1 z1]
 *  |   ...  |
 *  [xn yn zn] 
 * 
 * mat must be initialized. */
int sift3d_keypoint_store_to_mat_rm(const sift3d_keypoint_store *const kp,
                                    sift3d_mat_rm *const mat) {
    int i;

    const int num = kp->slab.num;

    // Resize mat
    mat->num_rows = num;
    mat->num_cols = IM_NDIMS;
    mat->type = SIFT3D_DOUBLE;
    if (resize_Mat_rm(mat))
        return SIFT3D_FAILURE;

    // Build the matrix
    for (i = 0; i < num; i++) {

        const sift3d_keypoint *const key = kp->buf + i;

        // Adjust the coordinates to the base octave
        const double coord_factor = ldexp(1.0, key->o);

        SIFT3D_MAT_RM_GET(mat, i, 0, double) = coord_factor * key->xd;
        SIFT3D_MAT_RM_GET(mat, i, 1, double) = coord_factor * key->yd;
        SIFT3D_MAT_RM_GET(mat, i, 2, double) = coord_factor * key->zd;
    }

    return SIFT3D_SUCCESS;
}

/* Convert SIFT3D descriptors to a matrix.
 *
 * Output format:
 *  [x y z el0 el1 ... el767]
 * Each row is a feature descriptor. [x y z] are the coordinates in image
 * space, and [el0 el1 ... el767 are the 768 dimensions of the descriptor.
 *
 * mat must be initialized prior to calling this function. mat will be resized.
 * The type of mat will be changed to float.
 */
int sift3d_descriptor_store_to_mat_rm(const sift3d_descriptor_store *const store,
                                      sift3d_mat_rm *const mat) {
    int i, j, a, p;

    const int num_rows = store->num;
    const int num_cols = IM_NDIMS + DESC_NUMEL;

    // Verify inputs
    if (num_rows < 1) {
        printf("SIFT3D_Descriptor_store_to_Mat_rm: invalid number of "
               "descriptors: %d \n", num_rows);
        return SIFT3D_FAILURE;
    }

    // Resize inputs
    mat->type = SIFT3D_FLOAT;
    mat->num_rows = num_rows;
    mat->num_cols = num_cols;
    if (resize_Mat_rm(mat))
        return SIFT3D_FAILURE;

    // Copy the data
    for (i = 0; i < num_rows; i++) {

        const sift3d_descriptor *const desc = store->buf + i;

        // Copy the coordinates
        SIFT3D_MAT_RM_GET(mat, i, 0, float) = (float) desc->xd;
        SIFT3D_MAT_RM_GET(mat, i, 1, float) = (float) desc->yd;
        SIFT3D_MAT_RM_GET(mat, i, 2, float) = (float) desc->zd;

        // Copy the feature vector
        for (j = 0; j < DESC_NUM_TOTAL_HIST; j++) {
            const sift3d_hist *const hist = desc->hists + j;
            HIST_LOOP_START(a, p)
                const int col = DESC_MAT_GET_COL(j, a, p);
            SIFT3D_MAT_RM_GET(mat, i, col, float) = 
                HIST_GET(hist, a, p);
            HIST_LOOP_END
                }
    }

    return SIFT3D_SUCCESS;
}

/* Write a sift3d_keypoint_store to a text file. The keypoints are stored in a matrix
 * (.csv, .csv.gz), where each keypoint is a row. The elements of each row are
 * as follows:
 *
 * x y z o s ori11 ori12 ... orinn
 *
 * x - the x-coordinate
 * y - the y-coordinate
 * z - the z-coordinate
 * o - the pyramid octave. To convert to image coordinates, multiply x,y,z by 
 *      pow(2, o)
 * s - the scale coordinate
 * ori(ij) - the ith row, jth column of the orientation matrix */
int sift3d_keypoint_store_save(const char *path, const sift3d_keypoint_store *const kp) {
    sift3d_mat_rm mat;
    int i, i_R, j_R;

    // sift3d_keypoint data format constants
    const int kp_str = 0; // column of x-coordinate
    const int kp_x   = 1; // column of x-coordinate
    const int kp_y   = 2; // column of y-coordinate
    const int kp_z   = 3; // column of z-coordinate
    const int kp_o   = 4; // column of octave index
    const int kp_s   = 5; // column of s-coordinate
    const int kp_ori = 6; // first column of the orientation matrix
    const int ori_numel = IM_NDIMS * IM_NDIMS; // Number of orientation 
    // elements
    const int num_rows = kp->slab.num;
    const int num_cols = kp_ori + ori_numel;

    // Initialize the matrix
    if (init_Mat_rm(&mat, num_rows, num_cols, SIFT3D_DOUBLE, 
                    SIFT3D_FALSE))
        return SIFT3D_FAILURE;
       
    // Write the keypoints 
    for (i = 0; i < num_rows; i++) {

        const sift3d_keypoint *const key = kp->buf + i;
        const sift3d_mat_rm *const R = &key->R;

        // Write strength
        SIFT3D_MAT_RM_GET(&mat, i, kp_str, double) = key->strength;

        // Write the coordinates 
        SIFT3D_MAT_RM_GET(&mat, i, kp_x, double) = key->xd;
        SIFT3D_MAT_RM_GET(&mat, i, kp_y, double) = key->yd; 
        SIFT3D_MAT_RM_GET(&mat, i, kp_z, double) = key->zd; 
        SIFT3D_MAT_RM_GET(&mat, i, kp_o, double) = key->o; 
        SIFT3D_MAT_RM_GET(&mat, i, kp_s, double) = key->sd;

        // Write the orientation matrix
        SIFT3D_MAT_RM_LOOP_START(R, i_R, j_R)

            const int kp_idx = kp_ori + 
            SIFT3D_MAT_RM_GET_IDX(R, i_R, j_R);
        
        SIFT3D_MAT_RM_GET(&mat, i, kp_idx, double) = 
            (double) SIFT3D_MAT_RM_GET(R, i_R, j_R, float);

        SIFT3D_MAT_RM_LOOP_END
            }

    // Write the matrix 
    if (write_Mat_rm(path, &mat))
        goto write_kp_quit;

    // Clean up
    cleanup_Mat_rm(&mat);

    return SIFT3D_SUCCESS;

write_kp_quit:
    cleanup_Mat_rm(&mat);
    return SIFT3D_FAILURE;
}

/* Write sift3d_detector descriptors to a text file.
 * See SIFT3D_Descriptor_store_to_sift3d_mat_rm for the file format. */
int sift3d_descriptor_store_save(const char *path,
                                 const sift3d_descriptor_store *const desc) {
    sift3d_mat_rm mat;

    // Initialize the matrix
    if (init_Mat_rm(&mat, 0, 0, SIFT3D_FLOAT, SIFT3D_FALSE))
        return SIFT3D_FAILURE;
     
    // Write the data into the matrix 
    if (sift3d_descriptor_store_to_mat_rm(desc, &mat))
        goto write_desc_quit;

    // Write the matrix to the file
    if (write_Mat_rm(path, &mat))
        goto write_desc_quit;

    // Clean up
    cleanup_Mat_rm(&mat);
    return SIFT3D_SUCCESS;

write_desc_quit:
    cleanup_Mat_rm(&mat);
    return SIFT3D_FAILURE;
}

static int keypoint_strength_cmp (const void *p1, const void *p2) {
    const sift3d_keypoint *kp1 = p1;
    const sift3d_keypoint *kp2 = p2;

    return (kp1->strength < kp2->strength) ? 1 : -1;
}

// Public API

sift3d_detector* sift3d_make_detector() {
    sift3d_detector *detector = malloc (sizeof (sift3d_detector));

    if (init_SIFT3D (detector)) {
        // FIXME: Is it safe to cleanup it? Hell I know
        goto failure;
    }

    return detector;

failure:
    sift3d_free_detector (detector);
    return NULL;
}

void sift3d_free_detector(sift3d_detector *detector) {
    cleanup_SIFT3D (detector);
    free (detector);
}

sift3d_keypoint_store* sift3d_make_keypoint_store() {
    sift3d_keypoint_store *store = malloc (sizeof (sift3d_keypoint_store));
    init_Keypoint_store (store);

    return store;
}

void sift3d_free_keypoint_store(sift3d_keypoint_store *store) {
    cleanup_Keypoint_store (store);
    free (store);
}

sift3d_descriptor_store* sift3d_make_descriptor_store() {
    sift3d_descriptor_store *store = malloc (sizeof (sift3d_descriptor_store));
    init_SIFT3D_Descriptor_store (store);

    return store;
}

void sift3d_free_descriptor_store(sift3d_descriptor_store *store) {
    cleanup_SIFT3D_Descriptor_store (store);
    free (store);
}

void sift3d_keypoint_store_sort_by_strength (sift3d_keypoint_store *const store, int limit) {
    qsort (store->buf, store->slab.num, sizeof (sift3d_keypoint), keypoint_strength_cmp);

    if (store->slab.num > limit && limit != 0) {
        resize_Keypoint_store (store, limit);
    }
}

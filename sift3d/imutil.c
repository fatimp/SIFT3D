/* -----------------------------------------------------------------------------
 * imutil.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Miscellaneous utility routines for image processing, linear algebra, and 
 * statistical regression. This library completely defines the Image,
 * Mat_rm, and Ransac types, among others, and stands apart from the other
 * source.
 * -----------------------------------------------------------------------------
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <zlib.h>
#include <imutil.h>
#include "immacros.h"
#include "imtypes_private.h"
#include "imutil_private.h"
#include "nifti.h"

/* Stringify a macro name */
#define STR(x) #x

/* Stringify the result of a macro expansion */
#define XSTR(x) STR(x)

/* zlib definitions */
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#include <fcntl.h>
#include <io.h>
#define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#define SET_BINARY_MODE(file)
#endif

/* Supported file extensions */
const char ext_analyze[] = "img";;
const char ext_gz[] = "gz";
const char ext_nii[] = "nii";

/* LAPACK declarations */
#ifdef SIFT3D_MEX
// Set the integer width to Matlab's defined width
#include <uchar.h>
#include "mex.h"
typedef mwSignedIndex fortran_int;
#ifdef _WINDOWS
// Remove underscores from FORTRAN functions
#define dsyevd_ dsyevd
#endif
#else
typedef int32_t fortran_int;
#endif
extern void dsyevd_(const char *, const char *, const fortran_int *, double *,
                    const fortran_int *, double *, double *,
                    const fortran_int *, fortran_int *, const fortran_int *,
                    fortran_int *);

/* As realloc, but frees the underlying pointer and returns NULL on error, or
 * if size is 0 and ptr is non-NULL. */
void *SIFT3D_safe_realloc(void *ptr, size_t size) {

    void *ret;

    // Call realloc and handle failures
    if (size == 0 || (ret = realloc(ptr, size)) == NULL) {
        if (ptr != NULL) {
            free(ptr);
        }
        return NULL;
    }

    return ret;

}

/* Initialize a triangle mesh for first use. This must be called before mesh
 * can be used in any other functions. */
void init_Mesh(sift3d_mesh * const mesh)
{
    mesh->tri = NULL;
    mesh->num = -1;
}

/* Release all memory associated with a triangle mesh. mesh cannot be reused
 * before it is reinitialized. */
void cleanup_Mesh(sift3d_mesh * const mesh)
{
    free(mesh->tri);
}

/* Set all elements to zero */
static int zero_Mat_rm(sift3d_mat_rm *const mat) {
    int i, j;

#define SET_ZERO(type)                                  \
    SIFT3D_MAT_RM_LOOP_START(mat, i, j)                 \
        SIFT3D_MAT_RM_GET(mat, i, j, type) = (type) 0;  \
    SIFT3D_MAT_RM_LOOP_END

    switch (mat->type) {
    case SIFT3D_DOUBLE:
        SET_ZERO(double);
        break;
    case SIFT3D_FLOAT:
        SET_ZERO(float);
        break;
    case SIFT3D_INT:
        SET_ZERO(int);
        break;
    default:
        return SIFT3D_FAILURE;
    }
#undef SET_ZERO

    return SIFT3D_SUCCESS;
}

/* Shortcut function to initalize a matrix.
 * 
 * Parameters:
 *      mat - The matrix to be initialized
 *      num_rows - The number of rows
 *      num_cols - The number of columns
 *      type - The data type 
 *      set_zero - If true, initializes the elements to zero.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int init_Mat_rm(sift3d_mat_rm *const mat, const int num_rows, const int num_cols,
                const sift3d_mat_type type, const int set_zero) {
    mat->type = type;
    mat->num_rows = num_rows;
    mat->num_cols = num_cols;
    mat->u.data_double = NULL;
    mat->size = 0;
    mat->static_mem = SIFT3D_FALSE;

    if (resize_Mat_rm(mat))
        return SIFT3D_FAILURE;
        
    if (set_zero && zero_Mat_rm(mat))
        return SIFT3D_FAILURE;
    
    return SIFT3D_SUCCESS;
}

/* As init_Mat_rm, but aliases data memory with pointer p. The flag 
 * mat->static_mem is set, and the matrix does not need to be freed with 
 * cleanup_Mat_rm. But, an error will be thrown if the user attempts to resize
 * the memory. That is, resize_sift3d_mat_rm will only return success if the size of 
 * the matrix does not change. */ 
int init_Mat_rm_p(sift3d_mat_rm *const mat, const void *const p, const int num_rows,
                  const int num_cols, const sift3d_mat_type type,
                  const int set_zero) {
    // Perform normal initialization
    if (init_Mat_rm(mat, num_rows, num_cols, type, set_zero))
        return SIFT3D_FAILURE;

    // Clean up any existing memory
    cleanup_Mat_rm(mat);

    // Alias with provided memory and set the static flag
    mat->u.data_double = (double *) p;
    mat->static_mem = SIFT3D_TRUE;

    // Optionally set to zero 
    if (set_zero && zero_Mat_rm(mat))
        return SIFT3D_FAILURE;

    return SIFT3D_SUCCESS;
}

/* Prints the type of mat into the string str. */
void sprint_type_Mat_rm(const sift3d_mat_rm * const mat, char *const str)
{
    switch (mat->type) {
    case SIFT3D_DOUBLE:
        sprintf(str, "double");
        break;
    case SIFT3D_FLOAT:
        sprintf(str, "float");
        break;
    case SIFT3D_INT:
        sprintf(str, "int");
        break;
    default:
        sprintf(str, "<sprint_type_Mat_rm: unknown type>");
    }
}

/* Copies a matrix. dst will be resized. */
int copy_Mat_rm(const sift3d_mat_rm * const src, sift3d_mat_rm * const dst)
{

    // Resize dst
    dst->type = src->type;
    dst->num_rows = src->num_rows;
    dst->num_cols = src->num_cols;
    if (resize_Mat_rm(dst))
        return SIFT3D_FAILURE;

    // Copy the data (use memmove because of static mode)
    memmove(dst->u.data_double, src->u.data_double, src->size);

    return SIFT3D_SUCCESS;
}

/* Re-sizes a matrix. The following fields
 * must already be initialized:
 * -num_rows
 * -num_cols
 * -type
 * -u.data_* (NULL for first use, non-null for resize)
 *
 * The following fields will be modified:
 * -size
 * -u.data_* (Change is not guaranteed)
 * 
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int resize_Mat_rm(sift3d_mat_rm *const mat) {

    size_t type_size, total_size;

    const int num_rows = mat->num_rows;
    const int num_cols = mat->num_cols;
    double **const data = &mat->u.data_double;
    const size_t numel = (size_t)num_rows * (size_t)num_cols;
    const sift3d_mat_type type = mat->type;

    // Get the size of the underyling datatype
    switch (type) {
    case SIFT3D_DOUBLE:
        type_size = sizeof(double);
        break;
    case SIFT3D_FLOAT:
        type_size = sizeof(float);
        break;
    case SIFT3D_INT:
        type_size = sizeof(int);
        break;
    default:
        SIFT3D_ERR("resize_Mat_rm: unknown type! \n");
        return SIFT3D_FAILURE;
    }

    // Calculate the new size
    total_size = type_size * numel;

    // Do nothing if the size has not changed
    if (total_size == mat->size)
        return SIFT3D_SUCCESS;
    mat->size = total_size;

    // Check for static reallocation
    if (mat->static_mem) {
        SIFT3D_ERR("resize_Mat_rm: illegal re-allocation of static matrix \n");
        return SIFT3D_FAILURE;
    }

    // Reset if the new size is 0 
    if (total_size == 0) {
        cleanup_Mat_rm(mat);
        return init_Mat_rm(mat, num_rows, num_cols, type, SIFT3D_FALSE);
    }

    // Re-allocate the memory
    if ((*data = (double *) SIFT3D_safe_realloc(*data, total_size)) == NULL) {
        mat->size = 0;
        return SIFT3D_FAILURE;
    }

    return SIFT3D_SUCCESS;
}

/* De-allocate the memory for a sift3d_mat_rm struct, unless it was initialized in
 * static mode. */
void cleanup_Mat_rm(sift3d_mat_rm *mat) {

    if (mat->u.data_double == NULL)
        return;

    if (!mat->static_mem)
        free(mat->u.data_double);
}

/* Separate the file name component from its path */
static const char *get_file_name(const char *path) {

    const char *name;

    // Get the last file separator
    name = strrchr(path, SIFT3D_FILE_SEP);
    return name == NULL ? path : name;
}

/* Get the extension of a file name */
static const char *get_file_ext(const char *name)
{

    const char *dot;

    // Get the file name component
    name = get_file_name(name);

    // Get the last dot
    dot = strrchr(name, '.');

    return dot == NULL || dot == name ? "" : dot + 1;
}

/* Detect the format of the supplied file name. */
static sift3d_im_format im_get_format(const char *path) {
    const char *ext;

    // If not a directory, get the file extension
    ext = get_file_ext(path);

    // Check the known types
    if (!strcmp(ext, ext_analyze) || !strcmp(ext, ext_gz) ||
        !strcmp(ext, ext_nii))
        return NIFTI;

    // The type was not recognized
    return UNKNOWN;
}

/* Read an image from a file. The file extension must match one of the
 * supported formats.
 *
 * Supported formats:
 * - Analyze (.img, .img.gz)
 * - Directory of DICOM files
 * - NIFTI-1 (.nii, .nii.gz)
 *
 * Return values:
 * -SIFT3D_SUCCESS - sift3d_image successfully read
 * -SIFT3D_FILE_DOES_NOT_EXIST - The file does not exist
 * -SIFT3D_UNSUPPORTED_FILE_TYPE - The file type is not supported
 * -SIFT3D_WRAPPER_NOT_COMPILED - The file type is supported, but the wrapper 
 *      library was not compiled.
 * -SIFT3D_UNEVEN_SPACING - The image slices are unevenly spaced.
 * -SIFT3D_INCONSISTENT_AXES  - The image slices have inconsistent axes.
 * -SIFT3D_DUPLICATE_SLICES - There are multiple slices in the same location.
 * -SIFT3D_FAILURE - Other error
 */
int im_read(const char *path, sift3d_image *const im) {
    int ret;

    // Get the file format and write the file
    switch (im_get_format(path)) {
    case ANALYZE:
    case NIFTI:
        ret = read_nii(path, im);
        break;
    case FILE_ERROR:
        ret = SIFT3D_FAILURE;
        break;
    case UNKNOWN:
    default:
        SIFT3D_ERR("im_read: unrecognized file extension "
                   "from file %s \n", path);
        ret = SIFT3D_UNSUPPORTED_FILE_TYPE;
    }

    return ret;
}

/* Write an image to a file.
 * 
 * Supported formats:
 * -Directory of DICOM files
 * -NIFTI (.nii, .nii.gz)
 *
 * Return values:
 * -SIFT3D_SUCCESS - Successfully wrote the image
 * -SIFT3D_UNSUPPORTED_FILE_TYPE - Cannot write this file type
 * -SIFT3D_FAILURE - Other error
 */
int im_write(const char *path, const sift3d_image *const im) {
    // Get the file format 
    switch (im_get_format(path)) {
    case ANALYZE:
    case NIFTI:
        return write_nii(path, im);
    case UNKNOWN:
    default:
        // Otherwise, the file extension was not found
        SIFT3D_ERR("im_write: unrecognized file extension "
                   "from file %s \n", path);

        return SIFT3D_UNSUPPORTED_FILE_TYPE;
    }

    // Unreachable code
    return SIFT3D_FAILURE;
}

/* Write a matrix to a .csv or .csv.gz file. */
int write_Mat_rm(const char *path, const sift3d_mat_rm * const mat)
{

    FILE *file;
    gzFile gz;
    const char *ext;
    int i, j, compress;

    const char *mode = "w";

    // Get the file extension
    ext = get_file_ext(path);

    // Check if we need to compress the file
    compress = strcmp(ext, ext_gz) == 0;

    // Open the file
    if (compress) {
        if ((gz = gzopen(path, mode)) == Z_NULL)
            return SIFT3D_FAILURE;
    } else {
        if ((file = fopen(path, mode)) == NULL)
            return SIFT3D_FAILURE;
    }

#define WRITE_MAT(mat, format, type)                            \
    SIFT3D_MAT_RM_LOOP_START(mat, i, j)                         \
        const char delim = j < mat->num_cols - 1 ? ',' : '\n';  \
    if (compress) {                                             \
        gzprintf(gz, format, SIFT3D_MAT_RM_GET(mat, i, j,       \
                                               type));          \
        gzputc(gz, delim);                                      \
    } else {                                                    \
        fprintf(file, format, SIFT3D_MAT_RM_GET(mat, i, j,      \
                                                type));         \
        fputc(delim, file);                                     \
    }                                                           \
    SIFT3D_MAT_RM_LOOP_END

    // Write the matrix
    switch (mat->type) {
    case SIFT3D_DOUBLE:
        WRITE_MAT(mat, "%f", double); 
        break;
    case SIFT3D_FLOAT:
        WRITE_MAT(mat, "%f", float); 
        break;
    case SIFT3D_INT:
        WRITE_MAT(mat, "%d", int); 
        break;
    default:
        goto write_mat_quit;
    }
#undef WRITE_MAT

    // Check for errors and finish writing the matrix
    if (compress) {
        if (gzclose(gz) != Z_OK)
            goto write_mat_quit;
    } else {
        if (ferror(file))
            goto write_mat_quit;
        fclose(file);
    }

    return SIFT3D_SUCCESS;

write_mat_quit:
    if (compress) {
        gzclose(gz);
    } else {
        fclose(file);
    }
    return SIFT3D_FAILURE;
}

/* Zero an image. */
static void im_zero(sift3d_image * im) {
    int x, y, z, c;

    SIFT3D_IM_LOOP_START_C(im, x, y, z, c)
        SIFT3D_IM_GET_VOX(im, x, y, z, c) = 0.0f;
    SIFT3D_IM_LOOP_END_C
}

/* Shortcut to initialize an image for first-time use.
 * Allocates memory, and assumes the default stride. This
 * function calls init_im and initializes all values to 0. */
int init_im_with_dims(sift3d_image *const im, const int nx, const int ny, const int nz,
                      const int nc)
{

    init_im(im);
    im->nx = nx;
    im->ny = ny;
    im->nz = nz;
    im->nc = nc;

    im_default_stride(im);
    if (im_resize(im))
        return SIFT3D_FAILURE;

    im_zero(im);

    return SIFT3D_SUCCESS;
}

/* Calculate the strides of an image object in the default
 * manner. The following parameters must be initialized:
 * -nx
 * -ny
 * -nz
 * -nc
 * If a dimension is not used, its size should be set
 * to 1. */
void im_default_stride(sift3d_image *const im)
{

    size_t prod;
    int i;

    prod = (size_t) im->nc;
    SIFT3D_IM_GET_STRIDES(im)[0] = prod;

    for (i = 1; i < IM_NDIMS; i++) {
        prod *= SIFT3D_IM_GET_DIMS(im)[i - 1];
        SIFT3D_IM_GET_STRIDES(im)[i] = prod;
    }
}

/* Resize an image according to the current nx, ny,
 * and nz. Does not modify scale space information or
 * strides. Prior to calling this function, use init_im(im)
 * and initialize the following fields:
 * -nx
 * -ny
 * -nz
 * -nc
 * -xs (can be set by im_default_stride(im)) 
 * -ys (can be set by im_default_stride(im)) 
 * -zs (can be set by im_default_stride(im)) 
 *
 * All of this initialization can also be done with
 * init_im_with_dims(), which calls this function.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int im_resize(sift3d_image *const im)
{
    int i;

	//FIXME: This will not work for strange strides
    const size_t size = (size_t)im->nx * (size_t)im->ny * (size_t)im->nz * (size_t)im->nc;

    // Verify inputs
    for (i = 0; i < IM_NDIMS; i++) {

        const int dim = SIFT3D_IM_GET_DIMS(im)[i];

        if (dim > 0)
            continue;

        SIFT3D_ERR("im_resize: invalid dimension %d: %d \n", i,
                   dim);
        return SIFT3D_FAILURE;
    }
    if (im->nc < 1) {
        SIFT3D_ERR("im_resize: invalid number of channels: %d \n",
                   im->nc);
        return SIFT3D_FAILURE;
    }

    // Do nothing if the size has not changed
    if (im->size == size)
        return SIFT3D_SUCCESS;
    im->size = size;

    // Allocate new memory
    im->data = SIFT3D_safe_realloc(im->data, size * sizeof(float));

    return size != 0 && im->data == NULL ? SIFT3D_FAILURE : SIFT3D_SUCCESS;
}

/* Downsample an image by a factor of 2 in each dimension.
 * This function initializes dst with the proper 
 * dimensions, and allocates memory. */
int im_downsample_2x(const sift3d_image *const src, sift3d_image *const dst)
{

    int x, y, z, c;

    // Initialize dst
    dst->nx = (int)floor((double)src->nx / 2.0);
    dst->ny = (int)floor((double)src->ny / 2.0);
    dst->nz = (int)floor((double)src->nz / 2.0);
    dst->nc = src->nc;
    im_default_stride(dst);
    if (im_resize(dst))
        return SIFT3D_FAILURE;

    // Downsample
    SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)

        const int src_x = x << 1;
    const int src_y = y << 1;
    const int src_z = z << 1;

    SIFT3D_IM_GET_VOX(dst, x, y, z, c) =
        SIFT3D_IM_GET_VOX(src, src_x, src_y, src_z, c);
    SIFT3D_IM_LOOP_END_C 

        return SIFT3D_SUCCESS;
}

/* Copy an image's dimensions and stride into another. 
 * This function resizes dst.
 * 
 * @param src The source image.
 * @param dst The destination image.
 * @return Returns SIFT3D_SUCCESS or SIFT3D_FAILURE.
 */
static int im_copy_dims(const sift3d_image * const src, sift3d_image * dst) {
    if (src->data == NULL)
        return SIFT3D_FAILURE;

    dst->nx = src->nx;
    dst->ny = src->ny;
    dst->nz = src->nz;
    dst->xs = src->xs;
    dst->ys = src->ys;
    dst->zs = src->zs;
    dst->nc = src->nc;
    dst->ux = src->ux;
    dst->uy = src->uy;
    dst->uz = src->uz;

    return im_resize(dst);
}

/* Copy an image's data into another. This function
 * changes the dimensions and stride of dst,
 * and allocates memory. */
int im_copy_data(const sift3d_image * const src, sift3d_image * const dst)
{

    int x, y, z, c;

    // Return if src has no data 
    if (src->data == NULL)
        return SIFT3D_FAILURE;

    // Return if src and dst are the same
    if (dst->data == src->data)
        return SIFT3D_SUCCESS;

    // Resize dst
    if (im_copy_dims(src, dst))
        return SIFT3D_FAILURE;

    // Copy data
    SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)
        SIFT3D_IM_GET_VOX(dst, x, y, z, c) = 
        SIFT3D_IM_GET_VOX(src, x, y, z, c);
    SIFT3D_IM_LOOP_END_C

        return SIFT3D_SUCCESS;
}

/* Clean up memory for an sift3d_image */
void im_free(sift3d_image * im)
{
    if (im->data != NULL)
        free(im->data);
}

/* Find the maximum absolute value of an image */
static float im_max_abs(const sift3d_image *const im) {

    float max;
    int x, y, z, c;

    max = 0.0f;
    SIFT3D_IM_LOOP_START_C(im, x, y, z, c)

        const float samp = fabsf(SIFT3D_IM_GET_VOX(im, x, y, z, c));
    max = SIFT3D_MAX(max, samp);

    SIFT3D_IM_LOOP_END_C

        return max;
}

/* Scale an image to the [-1, 1] range, where
 * the largest absolute value is 1. */
void im_scale(const sift3d_image *const im)
{

    int x, y, z, c;

    // Find the maximum absolute value
    const float max = im_max_abs(im);
    if (max == 0.0f)
        return;

    // Divide by the max 
    SIFT3D_IM_LOOP_START_C(im, x, y, z, c)
        SIFT3D_IM_GET_VOX(im, x, y, z, c) /= max;
    SIFT3D_IM_LOOP_END_C
        }

/* Subtract src2 from src1, saving the result in
 * dst.
 * Resizes dst. 
 */
int im_subtract(sift3d_image * src1, sift3d_image * src2, sift3d_image * dst)
{

    int x, y, z, c;

    // Verify inputs
    if (src1->nx != src2->nx ||
        src1->ny != src2->ny ||
        src1->nz != src2->nz || src1->nc != src2->nc)
        return SIFT3D_FAILURE;

    // Resize the output image
    if (im_copy_dims(src1, dst))
        return SIFT3D_FAILURE;

    SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)
        SIFT3D_IM_GET_VOX(dst, x, y, z, c) =
        SIFT3D_IM_GET_VOX(src1, x, y, z, c) -
        SIFT3D_IM_GET_VOX(src2, x, y, z, c);
    SIFT3D_IM_LOOP_END_C return SIFT3D_SUCCESS;
}

/* Convolve_sep for general filters */
static int convolve_sep_gen(const sift3d_image * const src,
                            sift3d_image * const dst, const sift3d_sep_fir_filter * const f,
                            const int dim, const double unit)
{
    register int x, y, z, c, d;

    register const int half_width = f->width / 2;
    register const int nx = src->nx;
    register const int ny = src->ny;
    register const int nz = src->nz;
    register const float conv_eps = 0.1f;
    register const int dim_end = SIFT3D_IM_GET_DIMS(src)[dim] - 1;
    register const float unit_factor =  unit /
        SIFT3D_IM_GET_UNITS(src)[dim];
    register const int unit_half_width = 
        (int) ceilf(half_width * unit_factor);
    int start[] = {0, 0, 0};
    int end[] = {nx - 1, ny - 1, nz - 1};

    // Compute starting and ending points for the convolution dimension
    start[dim] += unit_half_width;
    end[dim] -= unit_half_width + 1;

    //TODO: Convert this to convolve_x, which only convolves in x,
    // then make a wrapper to restride, transpose, convolve x, and transpose 
    // back

    // Resize the output, with the default stride
    if (im_copy_dims(src, dst))
        return SIFT3D_FAILURE;
    im_default_stride(dst);
    if (im_resize(dst))
        return SIFT3D_FAILURE;

    // Initialize the output to zeros
    im_zero(dst);

#define SAMP_AND_ACC(src, dst, tap, coords, c)                          \
    {                                                                   \
        float frac;                                                     \
                                                                        \
        const int idx_lo[] = {(coords)[0], (coords)[1], (coords)[2]};   \
        int idx_hi[] = {idx_lo[0], idx_lo[1], idx_lo[2]};               \
                                                                        \
        /* Convert the physical coordinates to integer indices*/        \
        idx_hi[dim] += 1;                                               \
        frac = (coords)[dim] - (float) idx_lo[dim];                     \
                                                                        \
        /* Sample with linear interpolation */                          \
        SIFT3D_IM_GET_VOX(dst, x, y, z, c) += (tap) *                   \
            ((1.0f - frac) *                                            \
             SIFT3D_IM_GET_VOX(src, idx_lo[0], idx_lo[1], idx_lo[2], c) + \
             frac *                                                     \
             SIFT3D_IM_GET_VOX(src, idx_hi[0], idx_hi[1], idx_hi[2], c)); \
    }

    // First pass: process the interior
#pragma omp parallel for private(x) private(y) private(c)
    SIFT3D_IM_LOOP_LIMITED_START_C(dst, x, y, z, c, start[0], end[0], 
                                   start[1], end[1], start[2], end[2])

        float coords[] = { x, y, z };

    for (d = -half_width; d <= half_width; d++) {

        const float tap = f->kernel[d + half_width];
        const float step = d * unit_factor;

        // Adjust the sampling coordinates
        coords[dim] -= step;

        // Sample
        SAMP_AND_ACC(src, dst, tap, coords, c);

        // Reset the sampling coordinates
        coords[dim] += step;
    }

    SIFT3D_IM_LOOP_END_C

        // Second pass: process the boundaries
#pragma omp parallel for private(x) private(y) private(c)
        SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)

        const int i_coords[] = { x, y, z };

    // Skip pixels we have already processed
    if (i_coords[dim] >= start[dim] && i_coords[dim] <= end[dim]) 
        continue;

    // Process the boundary pixel
    for (d = -half_width; d <= half_width; d++) {

        float coords[] = { x, y, z };
        const float tap = f->kernel[d + half_width];
        const float step = d * unit_factor;

        // Adjust the sampling coordinates
        coords[dim] -= step;

        // Mirror coordinates
        if ((int) coords[dim] < 0) {
            coords[dim] = -coords[dim];
            assert((int) coords[dim] >= 0);
        } else if ((int) coords[dim] >= dim_end) {
            coords[dim] = 2.0f * dim_end - coords[dim] -    
                conv_eps;
            assert((int) coords[dim] < dim_end);
        }

        // Sample
        SAMP_AND_ACC(src, dst, tap, coords, c);
    }

    SIFT3D_IM_LOOP_END_C 

#undef SAMP_AND_ACC

        return SIFT3D_SUCCESS;
}

/* Convolve_sep for symmetric filters. */
static int convolve_sep_sym(const sift3d_image * const src, sift3d_image * const dst,
                            const sift3d_sep_fir_filter * const f, const int dim,
                            const double unit)
{

    // TODO: Symmetry-specific function
    return convolve_sep_gen(src, dst, f, dim, unit);
}

/* Horizontally convolves a separable filter with an image, 
 * on CPU. Currently only works in 3D.
 * 
 * This function chooses among the best variant of convolve_sep* based on
 * compilation options and filter parameters.
 * 
 * Parameters: 
 * src - input image (initialized)
 * dst - output image (initialized) 
 int x, y, z;
 * f - filter to be applied
 * dim - dimension in which to convolve
 * unit - the spacing of the filter coefficients
 */
static int convolve_sep(const sift3d_image * const src,
                        sift3d_image * const dst, const sift3d_sep_fir_filter * const f,
                        const int dim, const double unit) {
    return f->symmetric ? 
        convolve_sep_sym(src, dst, f, dim, unit) : 
        convolve_sep_gen(src, dst, f, dim, unit);
}

/* Permute the dimensions of an image.
 *
 * Arguments: 
 * src - input image (initialized)
 * dim1 - input permutation dimension (x = 0, y = 1, z = 2)
 * dim2 - output permutation dimension (x = 0, y = 1, z = 2)
 * dst - output image (initialized)
 * 
 * example:
 * im_permute(src, dst, 0, 1) -- permute x with y in src
 *                              and save to dst
 */
static int im_permute(const sift3d_image * const src,
                      const int dim1, const int dim2,
                      sift3d_image * const dst) {
    int x, y, z, c;

    // Verify inputs
    if (dim1 < 0 || dim2 < 0 || dim1 > 3 || dim2 > 3) {
        printf("im_permute: invalid dimensions: dim1 %d dim2 %d \n",
               dim1, dim2);
        return SIFT3D_FAILURE;
    }

    // Check for the trivial case
    if (dim1 == dim2) {
        return im_copy_data(src, dst);
    }

    // Permute the units
    memcpy(SIFT3D_IM_GET_UNITS(dst), SIFT3D_IM_GET_UNITS(src), 
           IM_NDIMS * sizeof(double));
    SIFT3D_IM_GET_UNITS(dst)[dim1] = SIFT3D_IM_GET_UNITS(src)[dim2];
    SIFT3D_IM_GET_UNITS(dst)[dim2] = SIFT3D_IM_GET_UNITS(src)[dim1];

    // Resize the output
    memcpy(SIFT3D_IM_GET_DIMS(dst), SIFT3D_IM_GET_DIMS(src), 
           IM_NDIMS * sizeof(int));
    SIFT3D_IM_GET_DIMS(dst)[dim1] = SIFT3D_IM_GET_DIMS(src)[dim2];
    SIFT3D_IM_GET_DIMS(dst)[dim2] = SIFT3D_IM_GET_DIMS(src)[dim1];
    dst->nc = src->nc;
    im_default_stride(dst);
    if (im_resize(dst))
        return SIFT3D_FAILURE;

    // Transpose the data
    SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)

        int src_coords[] = {x, y, z};
    int temp;

    // Permute the coordinates
    temp = src_coords[dim1];
    src_coords[dim1] = src_coords[dim2];
    src_coords[dim2] = temp;

    // Copy the datum
    SIFT3D_IM_GET_VOX(dst, x, y, z, c) =
        SIFT3D_IM_GET_VOX(src, src_coords[0], src_coords[1], src_coords[2], c);

    SIFT3D_IM_LOOP_END_C 

        return SIFT3D_SUCCESS;
}

/* Computes the eigendecomposition of a real symmetric matrix, 
 * A = Q * diag(L) * Q', where Q is a real orthogonal matrix and L is a real 
 * diagonal matrix.
 *
 * A must be an [nxn] matrix. Q is [nxm], where m is in the interval [1, n],
 * depending on the values of A. L is [nx1], where the first m elements are
 * sorted in ascending order. The remaining n - m elements are zero. 
 * 
 * If Q is NULL, the eigenvectors will not be computed.
 *
 * The eigendecomposition is computed by divide and conquer.
 * 
 * This function resizes all non-null outputs and sets their type to double.
 *
 * This function does not ensure that A is symmetric.
 *
 * All matrices must be initialized prior to calling this funciton.
 * All matrices must have type double.
 *
 * Note: This function computes all of the eigenvalues, to a high degree of 
 * accuracy. A faster implementation is possible if you do not need high
 * precision, or if you do not need all of the eigenvalues, or if you do not 
 * need eigenvalues outside of some interval. 
 */
int eigen_Mat_rm(sift3d_mat_rm * A, sift3d_mat_rm * Q, sift3d_mat_rm * L)
{

    sift3d_mat_rm A_trans;
    double *work;
    fortran_int *iwork;
    double lwork_ret;
    fortran_int info, lwork, liwork;

    const char jobz = Q == NULL ? 'N' : 'V';
    const char uplo = 'U';
    const fortran_int n = A->num_cols;
    const fortran_int lda = n;
    const fortran_int lwork_query = -1;
    const fortran_int liwork_query = -1;

    // Verify inputs
    if (A->num_rows != n) {
        puts("eigen_Mat_rm: A be square \n");
        return SIFT3D_FAILURE;
    }
    if (A->type != SIFT3D_DOUBLE) {
        puts("eigen_Mat_rm: A must have type double \n");
        return SIFT3D_FAILURE;
    }
    // Resize outputs
    L->num_rows = n;
    L->num_cols = 1;
    L->type = SIFT3D_DOUBLE;
    if (resize_Mat_rm(L))
        return SIFT3D_FAILURE;

    // Initialize intermediate matrices and buffers
    work = NULL;
    iwork = NULL;
    if (init_Mat_rm(&A_trans, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE))
        goto EIGEN_MAT_RM_QUIT;

    // Copy the input matrix (A = A')
    if (copy_Mat_rm(A, &A_trans))
        goto EIGEN_MAT_RM_QUIT;

    // Query for the workspace sizes
    dsyevd_(&jobz, &uplo, &n, A_trans.u.data_double, &lda, L->u.data_double,
            &lwork_ret, &lwork_query, &liwork, &liwork_query, &info);

    if ((int32_t) info) {
        printf
            ("eigen_Mat_rm: LAPACK dsyevd workspace query error code %d",
             info);
        goto EIGEN_MAT_RM_QUIT;
    }
    // Allocate work spaces 
    lwork = (fortran_int) lwork_ret;
    if ((work = (double *)malloc(lwork * sizeof(double))) == NULL ||
        (iwork =
         (fortran_int *) malloc(liwork * sizeof(fortran_int))) == NULL)
        goto EIGEN_MAT_RM_QUIT;

    // Compute the eigendecomposition
    dsyevd_(&jobz, &uplo, &n, A_trans.u.data_double, &lda, L->u.data_double,
            work, &lwork, iwork, &liwork, &info);

    if ((int32_t) info) {
        printf("eigen_Mat_rm: LAPACK dsyevd error code %d", (int) info);
        goto EIGEN_MAT_RM_QUIT;
    }
    // Optionally return the eigenvectors
    if (Q != NULL && transpose_Mat_rm(&A_trans, Q))
        goto EIGEN_MAT_RM_QUIT;

    free(work);
    free(iwork);
    cleanup_Mat_rm(&A_trans);
    return SIFT3D_SUCCESS;

EIGEN_MAT_RM_QUIT:
    if (work != NULL)
        free(work);
    if (iwork != NULL)
        free(iwork);
    cleanup_Mat_rm(&A_trans);
    return SIFT3D_FAILURE;
}

/* Tranposes a matrix. Resizes dst with the type of src. 
 * All matrices must be initialized prior to calling this function. */
int transpose_Mat_rm(const sift3d_mat_rm *const src, sift3d_mat_rm *const dst)
{

    int i, j;

    // Verify inputs
    if (src->num_rows < 1 || src->num_cols < 1)
        return SIFT3D_FAILURE;

    // Resize the output
    dst->type = src->type;
    dst->num_rows = src->num_cols;
    dst->num_cols = src->num_rows;
    if (resize_Mat_rm(dst))
        return SIFT3D_FAILURE;

#define TRANSPOSE_MAT_RM(type)                  \
    SIFT3D_MAT_RM_LOOP_START(src, i, j)         \
        SIFT3D_MAT_RM_GET(dst, j, i, type) =    \
        SIFT3D_MAT_RM_GET(src, i, j, type);     \
    SIFT3D_MAT_RM_LOOP_END

    // Transpose
    switch (src->type) {
    case SIFT3D_DOUBLE:
        TRANSPOSE_MAT_RM(double);
        break;
    case SIFT3D_FLOAT:
        TRANSPOSE_MAT_RM(float);
        break;
    case SIFT3D_INT:
        TRANSPOSE_MAT_RM(int);
        break;
    default:
#ifndef NDEBUG
        puts("transpose_Mat_rm: unknown type \n");
#endif
        return SIFT3D_FAILURE;
    }
#undef TRANSPOSE_MAT_RM

    return SIFT3D_SUCCESS;
}

/* Apply a separable filter in multiple dimensions. This function resamples the
 * input to have the same units as f, then resamples the output to the
 * original units.
 *
 * Parameters:
 *  -src: The input image.
 *  -dst: The filtered image.
 *  -f: The filter to apply.
 *  -unit: The physical units of the filter kernel. Use -1.0 for the default,
 *      which is the same units as src.
 *
 * Return: SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int apply_Sep_FIR_filter(const sift3d_image * const src, sift3d_image * const dst,
                         sift3d_sep_fir_filter * const f, const double unit)
{

    sift3d_image temp;
    sift3d_image *cur_src, *cur_dst;
    int i;

    const double unit_default = -1.0;

    // Verify inputs
    if (unit < 0 && unit != unit_default) {
        SIFT3D_ERR("apply_Sep_FIR_filter: invalid unit: %f, use "
                   "%f for default \n", unit, unit_default);
        return SIFT3D_FAILURE;
    }

    // Resize the output
    if (im_copy_dims(src, dst))
        return SIFT3D_FAILURE; 

    // Allocate temporary storage
    init_im(&temp);
    if (im_copy_data(src, &temp))
        goto apply_sep_f_quit;

#define SWAP_BUFFERS                            \
    if (cur_dst == &temp) {                     \
        cur_src = &temp;                        \
        cur_dst = dst;                          \
    } else {                                    \
        cur_src = dst;                          \
        cur_dst = &temp;                        \
    }

    // Apply in n dimensions
    cur_src = (sift3d_image *) src;
    cur_dst = &temp;
    for (i = 0; i < IM_NDIMS; i++) {

        // Check for default parameters
        const double unit_arg = unit == unit_default ?
            SIFT3D_IM_GET_UNITS(src)[i] : unit;

        // Transpose so that the filter dimension is x
        if (i != 0) {
            if (im_permute(cur_src, 0, i, cur_dst))
                goto apply_sep_f_quit;
            SWAP_BUFFERS
                }

        // Apply the filter
        convolve_sep(cur_src, cur_dst, f, 0, unit_arg);
        SWAP_BUFFERS

            // Transpose back
            if (i != 0) {
                if (im_permute(cur_src, 0, i, cur_dst))
                    goto apply_sep_f_quit;
                SWAP_BUFFERS
            }
    }

    // Swap back
    SWAP_BUFFERS;

#undef SWAP_BUFFERS

    // Copy result to dst, if necessary
    if (cur_dst != dst && im_copy_data(cur_dst, dst))
        goto apply_sep_f_quit;

    // Clean up
    im_free(&temp);
    return SIFT3D_SUCCESS;

apply_sep_f_quit:
    im_free(&temp);
    return SIFT3D_FAILURE;
}

/* Initialize a separable FIR filter struct with the given parameters. If OpenCL
 * support is enabled and initialized, this creates a program to apply it with
 * separable filters.  
 *
 * Note that the kernel data will be copied, so the user can free it without 
 * affecting f. */
static int init_Sep_FIR_filter(sift3d_sep_fir_filter *const f,
                               const int dim, const int width,
                               const float *const kernel, const int symmetric) {
    const size_t kernel_size = width * sizeof(float);

    // Save the data
    f->dim = dim;
    f->width = width;
    f->symmetric = symmetric;

    // Allocate the kernel memory
    if ((f->kernel = (float *) malloc(kernel_size)) == NULL) {
        SIFT3D_ERR("init_Sep_FIT_filter: out of memory! \n");
        return SIFT3D_FAILURE;
    }

    // Copy the kernel data
    memcpy(f->kernel, kernel, kernel_size);
    return SIFT3D_SUCCESS;
}

/* Free a Sep_FIR_Filter. */
static void cleanup_Sep_FIR_filter(sift3d_sep_fir_filter *const f) {
    if (f->kernel != NULL) {
        free(f->kernel);
        f->kernel = NULL;
    }
}

/* Initialize the values of im so that it can be used by the
 * resize function. Does not allocate memory. */
void init_im(sift3d_image *const im)
{
    im->data = NULL;

    im->ux = 1;
    im->uy = 1;
    im->uz = 1;

    im->size = 0;
    im->s = -1.0;
    memset(SIFT3D_IM_GET_DIMS(im), 0, IM_NDIMS * sizeof(int));
    memset(SIFT3D_IM_GET_STRIDES(im), 0, IM_NDIMS * sizeof(size_t));
}

/* Initialize a normalized Gaussian filter, of the given sigma.
 * If SIFT3D_GAUSS_WIDTH_FCTR is defined, use that value for
 * the ratio between the width of the filter and sigma. Otherwise,
 * use the default value 3.0 
 */
#ifndef SIFT3D_GAUSS_WIDTH_FCTR
#define SIFT3D_GAUSS_WIDTH_FCTR 3.0
#endif
static int init_Gauss_filter(sift3d_gauss_filter * const gauss,
                             const double sigma,
                             const int dim) {
    float *kernel;
    double x;
    float acc;
    int i;

    const int half_width = sigma > 0 ? 
        SIFT3D_MAX((int)ceil(sigma * SIFT3D_GAUSS_WIDTH_FCTR), 1) :
        1;
    const int width = 2 * half_width + 1;

    // Initialize intermediates 
    if ((kernel = (float *) malloc(width * sizeof(float))) == NULL)
        return SIFT3D_FAILURE;

    // Calculate coefficients
    acc = 0;
    for (i = 0; i < width; i++) {
        // distance away from center of filter
        x = (double)i - half_width;

        // (x / sigma)^2 = x*x / (sigma*sigma)
        x /= sigma + DBL_EPSILON;

        // exponentiate result
        kernel[i] = (float)exp(-0.5 * x * x);

        // sum of all kernel elements
        acc += kernel[i];
    }

    // normalize kernel to sum to 1
    for (i = 0; i < width; i++) {
        kernel[i] /= acc;
    }

    // Save the filter data 
    gauss->sigma = sigma;
    if (init_Sep_FIR_filter(&gauss->f, dim, width, kernel,
                            SIFT3D_TRUE))
        goto init_Gauss_filter_quit;

    // Clean up
    free(kernel);

    return SIFT3D_SUCCESS;

init_Gauss_filter_quit:
    free(kernel);
    return SIFT3D_FAILURE;
}

/* Initialize a Gaussian filter to go from scale s_cur to s_next. */
static int init_Gauss_incremental_filter(sift3d_gauss_filter * const gauss,
                                         const double s_cur,
                                         const double s_next,
                                         const int dim) {
    double sigma;

    if (s_cur > s_next) {
        SIFT3D_ERR("init_Gauss_incremental_filter: "
                   "s_cur (%f) > s_next (%f) \n", s_cur, s_next);
        return SIFT3D_FAILURE;
    }
    assert(dim > 0);

    // Compute filter width parameter (sigma)
    sigma = sqrt(s_next * s_next - s_cur * s_cur);

    // Initialize filter kernel
    if (init_Gauss_filter(gauss, sigma, dim))
        return SIFT3D_FAILURE;

    return SIFT3D_SUCCESS;
}

/* Free a sift3d_gauss_filter */
static void cleanup_Gauss_filter(sift3d_gauss_filter * gauss) {
    cleanup_Sep_FIR_filter(&gauss->f);
}

/* Initialize a GSS filters stuct. This must be called before gss can be
 * used in any other functions. */
void init_GSS_filters(sift3d_gss_filters * const gss)
{
    gss->num_filters = -1;
    gss->gauss_octave = NULL;
}

/* Create GSS filters to create the given scale-space 
 * pyramid. */
int make_gss(sift3d_gss_filters * const gss, const sift3d_pyramid * const pyr)
{
    sift3d_image *cur, *next;
    int o, s;

    const int dim = 3;

    const int num_filters = pyr->num_levels - 1;
    const int first_level = pyr->first_level;
    const int last_level = SIFT3D_PYR_LAST_LEVEL(pyr);

    // Verify inputs
    if (num_filters < 1) {
        SIFT3D_ERR("make_gss: pyr has only %d levels, must have "
                   "at least 2", pyr->num_levels);
        return SIFT3D_FAILURE;
    }

    // Free all previous data, if any
    cleanup_GSS_filters(gss);
    init_GSS_filters(gss);

    // Copy pyramid parameters
    gss->num_filters = num_filters;
    gss->first_level = first_level;

    // Allocate the filter array (num_filters cannot be zero)
    if ((gss->gauss_octave = (sift3d_gauss_filter *) 
         SIFT3D_safe_realloc(gss->gauss_octave, 
                             num_filters * sizeof(sift3d_gauss_filter))) == NULL)
        return SIFT3D_FAILURE;

    // Make the filter for the very first blur
    next = SIFT3D_PYR_IM_GET(pyr, pyr->first_octave, first_level);
    if (init_Gauss_incremental_filter(&gss->first_gauss, pyr->sigma_n,
                                      next->s, dim))
        return SIFT3D_FAILURE;

    // Make one octave of filters (num_levels - 1)
    o = pyr->first_octave;
    for (s = first_level; s < last_level; s++) {
        cur = SIFT3D_PYR_IM_GET(pyr, o, s);
        next = SIFT3D_PYR_IM_GET(pyr, o, s + 1);
        if (init_Gauss_incremental_filter(SIFT3D_GAUSS_GET(gss, s),
                                          cur->s, next->s, dim))
            return SIFT3D_FAILURE;
    }

    return SIFT3D_SUCCESS;
}

/* Free all memory associated with the GSS filters. gss cannot be reused
 * unless it is reinitialized. */
void cleanup_GSS_filters(sift3d_gss_filters * const gss)
{
    int i;

    const int num_filters = gss->num_filters;

    // We are done if gss has no filters
    if (num_filters < 1)
        return;

    // Free the first filter
    cleanup_Gauss_filter(&gss->first_gauss);

    // Free the octave filters
    for (i = 0; i < num_filters; i++) {
        sift3d_gauss_filter *const g = gss->gauss_octave + i;
        cleanup_Gauss_filter(g);
    }

    // Free the octave filter buffer
    free(gss->gauss_octave);
}

/* Initialize a sift3d_pyramid for use. Must be called before a sift3d_pyramid can be used
 * in any other functions. */
void init_Pyramid(sift3d_pyramid * const pyr)
{
    pyr->levels = NULL;
    pyr->first_level = 0;
    pyr->num_levels = pyr->num_kp_levels = 0;
    pyr->first_octave = 0;
    pyr->num_octaves = 0;
    pyr->sigma0 = pyr->sigma_n = 0.0;
}

/* Resize a scale-space pyramid according to the size of base image im.
 *
 * Parameters:
 *  -im: An image with the desired dimensions and units at octave 0
 *  -first_level: The index of the first pyramid level per octave
 *  -num_kp_levels: The number of levels per octave in which keypoints are 
 *      detected
 *  -num_levels: The total number of levels. Must be greater than or equal to
 *      num_kp_levels.
 *  -first_octave: The index of the first octave (0 is the base)
 *  -num_octaves: The total number of octaves 
 *  -sigma0: The scale parameter of level 0, octave 0
 *  -sigma_n: The nominal scale of the image im.
 *  -pyr: The sift3d_pyramid to be resized.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int resize_Pyramid(const sift3d_image *const im, const int first_level, 
                   const unsigned int num_kp_levels, const unsigned int num_levels,
                   const int first_octave, const unsigned int num_octaves, 
                   sift3d_pyramid *const pyr) {

    double units[IM_NDIMS];
    int dims[IM_NDIMS];
    double factor;
    int i, o, s;

    const double sigma0 = pyr->sigma0;
    const double sigma_n = pyr->sigma_n;
    const int old_num_total_levels = pyr->num_levels * pyr->num_octaves;
    const int num_total_levels = num_levels * num_octaves;

    // Verify inputs
    if (num_levels < num_kp_levels) {
        SIFT3D_ERR("resize_Pyramid: num_levels (%u) < "
                   "num_kp_levels (%d)", num_levels, num_kp_levels);
        return SIFT3D_FAILURE;
    }

    // Store the new parameters
    pyr->first_level = first_level;
    pyr->num_kp_levels = num_kp_levels;
    pyr->first_octave = first_octave;
    pyr->num_octaves = num_octaves;
    pyr->num_levels = num_levels;

    // Clean up old levels which are no longer needed 
    for (i = num_total_levels; i < old_num_total_levels; i++) {
        sift3d_image *const level = pyr->levels + i;
        im_free(level);
    }

    // Resize the outer array
    if (num_total_levels != 0 && 
        ((pyr->levels =
          SIFT3D_safe_realloc(pyr->levels,
                              num_total_levels * sizeof(sift3d_image))) == NULL))
        return SIFT3D_FAILURE;

    // We have nothing more to do if there are no levels
    if (num_total_levels == 0)
        return SIFT3D_SUCCESS;

    // Initalize new levels
    for (i = old_num_total_levels; i < num_total_levels; i++) {
        sift3d_image *const level = pyr->levels + i;
        init_im(level);
    }

    // We have nothing more to do if the image is empty
    if (im->data == NULL)
        return SIFT3D_SUCCESS;

    // Calculate base image dimensions and units
    factor = pow(2.0, -first_octave);
    for (i = 0; i < IM_NDIMS; i++) {
        dims[i] = (int) ((double) SIFT3D_IM_GET_DIMS(im)[i] * factor);
        units[i] = SIFT3D_IM_GET_UNITS(im)[i] * factor;
    }

    // Initialize each level separately
    SIFT3D_PYR_LOOP_START(pyr, o, s)
        // Initialize sift3d_image fields
        sift3d_image *const level = SIFT3D_PYR_IM_GET(pyr, o, s);
    memcpy(SIFT3D_IM_GET_DIMS(level), dims, 
           IM_NDIMS * sizeof(int));
    memcpy(SIFT3D_IM_GET_UNITS(level), units, 
           IM_NDIMS * sizeof(double));
    level->nc = im->nc;
    im_default_stride(level);

    // Re-size data memory
    if (im_resize(level))
        return SIFT3D_FAILURE;

    SIFT3D_PYR_LOOP_SCALE_END

        // Adjust dimensions and recalculate image size
        for (i = 0; i < IM_NDIMS; i++) {
            dims[i] /= 2;
            units[i] *= 2;
        }

    SIFT3D_PYR_LOOP_OCTAVE_END 

        // Set the scales for the new levels
        return set_scales_Pyramid(pyr->sigma0, pyr->sigma_n, pyr);
}

/* Set the scale-space parameters on a sift3d_pyramid struct. Operates on all levels
 * of the pyramid. This function is called automatically by resize_Pyramid.
 *
 * Parameters:
 *  -sigma0: The scale parameter of level 0, octave 0
 *  -sigma_n: The nominal scale parameter of images being transfomed into
 *      this pyramid struct. 
 *  -Pyr: The sift3d_pyramid to be modified. */
int set_scales_Pyramid(const double sigma0, const double sigma_n, 
                       sift3d_pyramid *const pyr) {

    int o, s;

    const int num_kp_levels = pyr->num_kp_levels;
    const sift3d_image *const first_level = 
        SIFT3D_PYR_IM_GET(pyr, pyr->first_octave, pyr->first_level);

    // Compute the scales of each level
    SIFT3D_PYR_LOOP_START(pyr, o, s)

        // Compute the scale 
        sift3d_image *const level = SIFT3D_PYR_IM_GET(pyr, o, s);
    const double scale = 
        sigma0 * pow(2.0, o + (double) s / num_kp_levels);

    // Verify that sigma_n is not too large
    if (o == pyr->first_octave && s == pyr->first_level && 
        scale < sigma_n) {
        SIFT3D_ERR("set_scales_Pyramid: sigma_n too large "
                   "for these settings. Max allowed: %f \n", 
                   scale - DBL_EPSILON);
        return SIFT3D_FAILURE;
    }

    // Save the scale
    level->s = scale;
    SIFT3D_PYR_LOOP_END

        // Store the parameters
        pyr->sigma0 = sigma0;
    pyr->sigma_n = sigma_n;

    return SIFT3D_SUCCESS;
}

/* Release all memory associated with a Pyramid. pyr cannot be used again,
 * unless it is reinitialized. */
void cleanup_Pyramid(sift3d_pyramid * const pyr)
{

    int o, s;

    // We are done if there are no levels
    if (pyr->levels == NULL)
        return;

    // Free the levels
    SIFT3D_PYR_LOOP_START(pyr, o, s)
        sift3d_image *const level = SIFT3D_PYR_IM_GET(pyr, o, s);
    im_free(level);
    SIFT3D_PYR_LOOP_END

        // Free the pyramid level buffer
        free(pyr->levels);
}

/* Initialize a sift3d_slab for first use */
void init_Slab(sift3d_slab *const slab) {
    slab->buf_size = slab->num = 0;
    slab->buf = NULL;
}

/* Free all memory associated with a slab. sift3d_slab cannot be re-used after 
 * calling this function, unless re-initialized. */
void cleanup_Slab(sift3d_slab * const slab)
{
    if (slab->buf != NULL)
        free(slab->buf);
}


// Public API

void sift3d_free_image (sift3d_image *image) {
    im_free (image);
    free (image);
}

sift3d_image* sift3d_make_image (const int nx, const int ny, const int nz, const int nc) {
    sift3d_image *image = malloc (sizeof (sift3d_image));
    if (init_im_with_dims (image, nx, ny, nz, nc)) {
        goto failure;
    }

    return image;

failure:
    sift3d_free_image (image);
    return NULL;
}

sift3d_image* sift3d_read_image(const char *path) {
    sift3d_image *image = malloc (sizeof (sift3d_image));
    init_im (image);

    if (im_read (path, image)) {
        goto failure;
    }

    return image;

failure:
    sift3d_free_image (image);
    return NULL;
}

float* sift3d_image_data(const sift3d_image *image) {
    return image->data;
}

sift3d_mat_rm* sift3d_make_mat_rm() {
    sift3d_mat_rm *matrix = malloc (sizeof (sift3d_mat_rm));
    if (init_Mat_rm (matrix, 0, 0, SIFT3D_FLOAT, 0)) {
        goto failure;
    }

    return matrix;

failure:
    sift3d_free_mat_rm (matrix);
    return NULL;
}

void sift3d_free_mat_rm (sift3d_mat_rm *matrix) {
    cleanup_Mat_rm (matrix);
    free (matrix);
}

void* sift3d_mat_rm_data (sift3d_mat_rm *matrix) {
    return matrix->u.data_float;
}

void sift3d_mat_rm_dimensions (const sift3d_mat_rm *matrix, int *num_cols, int *num_rows) {
    if (num_cols != NULL) {
        *num_cols = matrix->num_cols;
    }

    if (num_rows != NULL) {
        *num_rows = matrix->num_rows;
    }
}

sift3d_mat_type sift3d_mat_rm_type (const sift3d_mat_rm *matrix) {
    return matrix->type;
}

# SIFT3D        {#mainpage}

## Introduction

SIFT3D is an analogue of the scale-invariant feature transform (SIFT) for
three-dimensional images. It leverages volumetric data and real-world units to
detect keypoints and extract a robust description of their content.

The original work is by Blaine Rister et al. and this is a refactored and
bugfixed version.

## Dependencies

To build SIFT3D yo'll need:

* LAPACK
* cmake (build-only dependency)
* nifticlib (optional, for examples)

## Typical workflow

This is a typical workflow. All error checks are omitted for clarity, but you
should check error codes returned by SIFT3D:

```c
#include <sift3d/imtypes.h>
#include <sift3d/imutils.h>
#include <sift3d/sift.h>
#include <string.h>

// The input must be in column-major format, which is unnatural for C,
// but the author of the original code chose this order.

void get_descriptors (float *array, int nx, int ny, int nz) {
    // At first we create an image object
    int nc = 1;
    sift3d_image *image = sift3d_make_image (nx, ny, nz, nc);
    float *data = sift3d_image_data (image);
    // Axial strides are guaranteed to be nc, nc * nx, nc * nx * ny
    // and the total number of elements is nc * nx * ny * nz
    memcpy (data, array, nx * ny * nz * sizeof (float));

    // Create a keypoint detector and a storage for keypoints and descriptors
    sift3d_detector *detector = sift3d_make_detector();
    sift3d_detector *kps = sift3d_make_keypoint_store();
    sift3d_detector *descs = sift3d_make_descriptor_store();

    // Detect keypoints
    sift3d_detect_keypoints (detector, image, kps);
    // Extract descriptors
    sift3d_extract_descriptors (detector, kps, descs);

    // Create a matrix
    sift3d_mat_rm *mat = sift3d_make_mat_rm();

    // Copy descriptors to a matrix
    sift3d_descriptor_store_to_mat_rm (descs, mat);

    // Get descriptors
    void *desc_data = sift3d_mat_rm_data (mat);
    // Get dimensions
    int cols, rows;
    sift3d_mat_rm_dimensions (mat, &cols, &rows);
    // Get element type, should be SIFT3D_FLOAT, desc_data can be cast to float*
    sift3d_mat_type type = sift3d_mat_rm_type (mat);

    // Do cleanup
    sift3d_free_image (image);
    sift3d_free_detector (detector);
    sift3d_free_descriptor_store (descs);
    sift3d_free_keypoint_store (kps);

    // Do what you want with the descriptors
    ;
    // Free the matrix
    sift3d_free_mat_rm (mat);
}
```

/* -----------------------------------------------------------------------------
 * kpSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains the CLI to extract SIFT3D keypoints and descriptors from
 * a single image. 
 * -----------------------------------------------------------------------------
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <imutil.h>
#include <sift.h>

/* Options */
#define KEYS 'a'
#define DESC 'b'
#define HELP 1

/* Help message */
const char help_msg[] = 
    "Usage: kpSift3D [image.nii] \n"
    "\n"
    "Detects SIFT3D keypoints and extracts their descriptors from an "
    "image.\n" 
    "\n"
    "Example: \n"
    " kpSift3D --keys keys.csv --desc desc.csv image.nii \n"
    "\n"
    "Output options: \n"
    " --keys [filename] \n"
    "       Specifies the output file name for the keypoints. \n"
    "       Supported file formats: .csv, .csv.gz \n"
    " --desc [filename] \n"
    "       Specifies the output file name for the descriptors. \n"
    "       Supported file formats: .csv, .csv.gz \n"
    "At least one of the output options must be specified. \n";

/* CLI for 3D SIFT */
int main(int argc, char *argv[]) {
    sift3d_image *image           = NULL;
    sift3d_detector *sift3d       = NULL;
    sift3d_keypoint_store *kp     = NULL;
    sift3d_descriptor_store *desc = NULL;
    int retcode                   = 0;
    char *im_path, *keys_path, *desc_path;
    int c, num_args;

    const struct option longopts[] = {
        {"keys", required_argument, NULL, KEYS},
        {"desc", required_argument, NULL, DESC},
        {"help", no_argument, NULL, HELP},
        {0, 0, 0, 0}
    };

    // Parse the kpSift3d options
    opterr = 1;
    keys_path = desc_path = NULL;
    while ((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch (c) {
        case KEYS:
            keys_path = optarg;
            break;
        case DESC:
            desc_path = optarg;
            break;
        case HELP:
            puts(help_msg);
            return 0;
        case '?':
        default:
            return 1;
        }
    }

    // Ensure we have at least one output
    if (keys_path == NULL && desc_path == NULL) {
        fprintf (stderr, "No outputs specified.\n");
        return 1;
    }

    // Parse the required arguments
    num_args = argc - optind;
    if (num_args < 1) {
        fprintf (stderr, "Not enough arguments.\n");
        return 1;
    } else if (num_args > 1) {
        fprintf (stderr, "Too many arguments.\n");
        return 1;
    }
    im_path = argv[optind];

    // Initialize the SIFT data 
    sift3d = sift3d_make_detector ();
    if (sift3d == NULL) {
        fprintf (stderr, "Failed to initialize SIFT data.\n");
        retcode = 1;
        goto done;
    }

    // Initialize data 
    kp = sift3d_make_keypoint_store();
    desc = sift3d_make_descriptor_store();

    // Read the image
    image = sift3d_read_image(im_path);
    if (image == NULL) {
        fprintf (stderr, "Could not read image.\n");
        retcode = 1;
        goto done;
    }

    // Extract keypoints
    if (sift3d_detect_keypoints(sift3d, image, kp)) {
        fprintf (stderr, "Failed to detect keypoints.\n");
        retcode = 1;
        goto done;
    }

    sift3d_keypoint_store_sort_by_strength (kp, 100);

    // Optionally write the keypoints 
    if (keys_path != NULL && sift3d_keypoint_store_save (keys_path, kp)) {
        fprintf(stderr, "Failed to write the keypoints to %s\n", keys_path);
        retcode = 1;
        goto done;
    }

    // Optionally extract descriptors
    if (desc_path != NULL) {

        // Extract descriptors
        if (sift3d_extract_descriptors(sift3d, kp, desc)) {
            fprintf (stderr, "Failed to extract descriptors.\n");
            retcode = 1;
            goto done;
        }

        // Write the descriptors
        if (sift3d_descriptor_store_save (desc_path, desc)) {
            fprintf(stderr, "Failed to write the descriptors to %s\n", desc_path);
            retcode = 1;
            goto done;
        }
    }

done:
    if (desc != NULL) {
        sift3d_free_descriptor_store (desc);
    }

    if (kp != NULL) {
        sift3d_free_keypoint_store (kp);
    }

    if (sift3d != NULL) {
        sift3d_free_detector (sift3d);
    }

    if (image != NULL) {
        sift3d_free_image (image);
    }

    return retcode;
}

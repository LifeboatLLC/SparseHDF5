/*
 * frame-writer-str.c
 *
 * This program creates an HDF5 file containing a 2D structured chunked dataset.
 *
 * Command-line options (all keyword=value):
 *   X=INT          - dataset dimension X
 *   Y=INT          - dataset dimension Y
 *   XC=INT         - chunk dimension in X
 *   YC=INT         - chunk dimension in Y
 *   m=FLOAT           - percentage of non-zero elements per chunk (0–100)
 *   compress=0|1      - whether to use gzip compression
 *   data=random|const - generate random or constant integer data
 *   pattern=random|contiguous - multiple random 1x1 hyperslabs or one contiguous
 *   outfile=STRING    - REQUIRED: output HDF5 file name
 *
 * For each chunk, the program:
 *   1. Defines a 2D hyperslab corresponding to that chunk.
 *   2. Generates data according to user options.
 *   3. Writes either:
 *        - one contiguous hyperslab, or
 *        - multiple small hyperslabs at random positions.
 *
 * Example:
 *   ./frame-writer-str X=100 Y=80 XC=20 YC=20 m=20 \
 *       compress=1 data=random pattern=random outfile=test.h5
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "hdf5.h"

#define DATASET_NAME "dset-str"

static void print_usage(const char *progname) {
    printf("Usage:\n");
    printf("  %s X=INT Y=INT XC=INT YC=INT m=FLOAT compress=0|1\n", progname);
    printf("     data=random|const pattern=random|contiguous outfile=FILENAME\n\n");
    printf("Example:\n");
    printf("  %s X=100 Y=80 XC=20 YC=20 m=20 compress=1 \\\n", progname);
    printf("     data=random pattern=contiguous outfile=test.h5\n");
}

static void parse_args(int argc, char *argv[],
                       hsize_t *xdim, hsize_t *ydim,
                       hsize_t *xchunk, hsize_t *ychunk,
                       double *mfrac, int *compress,
                       int *random_data, int *random_pattern,
                       char *outfile)
{
    int have_outfile = 0;

    *xdim = *ydim = 100;
    *xchunk = *ychunk = 10;
    *mfrac = 0.1;
    *compress = 0;
    *random_data = 1;
    *random_pattern = 1;
    outfile[0] = '\0';

    for (int i = 1; i < argc; i++) {
        char *arg = argv[i];
        if (strncmp(arg, "X=", 2) == 0)
            *xdim = atol(arg + 2);
        else if (strncmp(arg, "Y=", 2) == 0)
            *ydim = atol(arg + 2);
        else if (strncmp(arg, "XC=", 3) == 0)
            *xchunk = atol(arg + 3);
        else if (strncmp(arg, "YC=", 3) == 0)
            *ychunk = atol(arg + 3);
        else if (strncmp(arg, "m=", 2) == 0)
            *mfrac = atof(arg + 2) / 100.0;
        else if (strncmp(arg, "compress=", 9) == 0)
            *compress = atoi(arg + 9);
        else if (strncmp(arg, "data=", 5) == 0)
            *random_data = (strcmp(arg + 5, "random") == 0);
        else if (strncmp(arg, "pattern=", 8) == 0)
            *random_pattern = (strcmp(arg + 8, "random") == 0);
        else if (strncmp(arg, "outfile=", 8) == 0) {
            strcpy(outfile, arg + 8);
            have_outfile = 1;
        }
    }

    if (!have_outfile || outfile[0] == '\0') {
        fprintf(stderr, "Error: 'outfile' parameter is required.\n\n");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    hsize_t xdim, ydim, xchunk, ychunk;
    double mfrac;
    int compress, random_data, random_pattern;
    char outfile[256];
    unsigned int level        = 6;
    unsigned int cd_values[1] = {level};
    size_t       cd_nelmts    = 1;

    parse_args(argc, argv, &xdim, &ydim, &xchunk, &ychunk,
               &mfrac, &compress, &random_data, &random_pattern, outfile);

    srand((unsigned)time(NULL));

    hid_t file_id, space_id, dset_id, dcpl_id;
    hsize_t dims[2] = {xdim, ydim};
    hsize_t chunk_dims[2] = {xchunk, ychunk};
    hsize_t mem_dims[1];

    file_id = H5Fcreate(outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    space_id = H5Screate_simple(2, dims, NULL);
    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_struct_chunk(dcpl_id, 2, chunk_dims, H5D_SPARSE_CHUNK);
    if (compress) {
        H5Pset_filter2(dcpl_id, H5_SECTION_SELECTION, H5Z_FILTER_DEFLATE, H5Z_FLAG_OPTIONAL, cd_nelmts,
                           cd_values);
        H5Pset_filter2(dcpl_id, H5_SECTION_FIXED, H5Z_FILTER_DEFLATE, H5Z_FLAG_OPTIONAL, cd_nelmts,
                           cd_values);
    }

    dset_id = H5Dcreate(file_id, DATASET_NAME, H5T_NATIVE_INT, space_id,
                        H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    H5Sclose(space_id);
    H5Pclose(dcpl_id);


    hsize_t n_chunks_x = xdim / xchunk;
    hsize_t n_chunks_y = ydim / ychunk;

    for (hsize_t cx = 0; cx < n_chunks_x; cx++) {
        for (hsize_t cy = 0; cy < n_chunks_y; cy++) {
            hsize_t offset[2] = {cx * xchunk, cy * ychunk};
            hsize_t this_rows = xchunk;
            hsize_t this_cols = ychunk;

            hid_t filespace = H5Dget_space(dset_id);

            if (random_pattern) {
                // --- Random scattered 1x1 hyperslabs combined into one selection ---
                hsize_t num_points = (hsize_t)(this_rows * this_cols * mfrac);
                if (num_points == 0) num_points = 1;
                int *chunk_buf = (int *)malloc(num_points * sizeof(int));
                if (!chunk_buf) {
                    fprintf(stderr, "Memory allocation failed\n");
                    return 1;
                }

                // Create one combined selection and data buffer
                for (hsize_t n = 0; n < num_points; n++) {
                    hsize_t start[2] = {rand() % this_rows, rand() % this_cols};
                    //printf("Offsets within chunk %lld %lld \n", start[0], start[1]);
                    hsize_t slab[2] = {1, 1};
                    H5Sselect_hyperslab(filespace,
                                        (n == 0) ? H5S_SELECT_SET : H5S_SELECT_OR,
                                        (hsize_t[2]){offset[0] + start[0], offset[1] + start[1]},
                                        NULL, slab, NULL);

                    
                    /* Sometimes the same location is generated, to avoid it check if we
                       have the same number of selected points as iterations, if not, rerun */ 
                    while ((n+1)!= H5Sget_select_npoints(filespace)) {
                        hsize_t start[2] = {rand() % this_rows, rand() % this_cols};
                        H5Sselect_hyperslab(filespace, (n == 0) ? H5S_SELECT_SET : H5S_SELECT_OR,
                                        (hsize_t[2]){offset[0] + start[0], offset[1] + start[1]},
                                        NULL, slab, NULL);
                    } 
                    // Fill local buffer for visualization / test (optional)
                    chunk_buf[n] = random_data ? rand() % 100 : (int)(cx * n_chunks_y + cy + 1);
                }
                

                // Define full chunk memory space
                mem_dims[0] = num_points;
                hid_t memspace = H5Screate_simple(1, mem_dims, NULL);
                //printf ("Number of selected points %lld \n", num_points);
                hssize_t file_points = H5Sget_select_npoints(filespace);
                hssize_t mem_points = H5Sget_select_npoints(memspace);

                //printf ("Number of selected memory points %lld \n", mem_points);
                //printf ("Number of selected file points %lld \n", file_points);
                H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, chunk_buf);
                H5Sclose(memspace);
                H5Sselect_none(filespace);
                free(chunk_buf);
            } else {
                // --- One contiguous hyperslab per chunk ---
                //printf ("Percent %lf \n", mfrac);
                hsize_t slab_rows = (hsize_t)(this_rows * sqrt(mfrac));
                hsize_t slab_cols = (hsize_t)(this_cols * sqrt(mfrac));
                if (slab_rows == 0) slab_rows = 1;
                if (slab_cols == 0) slab_cols = 1;

                hsize_t start[2] = {rand() % (this_rows - slab_rows + 1),
                                    rand() % (this_cols - slab_cols + 1)};
                hsize_t slab[2] = {slab_rows, slab_cols};
                hsize_t total_elems = slab_rows * slab_cols;
                
                int *chunk_buf = (int *)malloc(total_elems * sizeof(int));
                if (!chunk_buf) {
                    fprintf(stderr, "Memory allocation failed\n");
                    return 1;
                }

                for (hsize_t i = 0; i < total_elems; i++)
                    chunk_buf[i] = random_data ? rand() % 100 : (int)(cx * n_chunks_y + cy + 1);

                hid_t memspace = H5Screate_simple(2, slab, NULL);
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET,
                                    (hsize_t[2]){offset[0] + start[0], offset[1] + start[1]},
                                    NULL, slab, NULL);
                H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, chunk_buf);
                H5Sclose(memspace);
                H5Sselect_none(filespace);
                free(chunk_buf);
            }

            H5Sclose(filespace);
        }
    }

    H5Dclose(dset_id);
    H5Fclose(file_id);

    printf("✅ HDF5 file '%s' created successfully.\n", outfile);
    return 0;
}


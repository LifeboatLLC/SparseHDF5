#include "hdf5.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FILE "sparse_write_and_read_chunk_test_minimal_hslab.h5"

#define SPARSE_DSET  "sparse_dset"

#define SUCCEED 1
#define FAIL -1

int main(void) {
    int ret_val = EXIT_SUCCESS;
    // char  filename[FILENAME_BUF_SIZE]; /* File name */
    hid_t   fid          = H5I_INVALID_HID;
    hid_t   sid          = H5I_INVALID_HID;
    hid_t   sid1         = H5I_INVALID_HID;
    hid_t   sid2         = H5I_INVALID_HID;
    hid_t   dcpl         = H5I_INVALID_HID;
    hid_t   did          = H5I_INVALID_HID;
    hsize_t dim[1]       = {5}; /* 1-d dataspace */
    hsize_t chunk_dim[1] = {5};  /* Chunk size */
    int     wbuf[5];            /* Write buffer (set to 5 to only have one chunk in the write) */ 
    herr_t  ret;
    int npoints;
    int rbuf[3];

    // Hyperslab selection components
    hid_t   fdataspace_id, fmemspace_id;
    hsize_t foffset[1] = {1};
    hsize_t fblock[1] = {1};
    hsize_t fcount[1] = {3};
    hsize_t fdims[1] = {3};
    hsize_t stride[1] = {1};
    int fdata[3];


    printf("APIs for direct chunk I/O on structured chunks\n");

    // SKIPPED();
    // return 0;

    /* Create a file */
    // h5_fixname(FILENAME[2], fapl, filename, sizeof filename);
    if ((fid = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        ret_val = EXIT_FAILURE;

    /* Create dataspace */
    if ((sid = H5Screate_simple(1, dim, NULL)) < 0)
        ret_val = EXIT_FAILURE;

    /* Create property list for compact dataset creation */
    if ((dcpl = H5Pcreate(H5P_DATASET_CREATE)) < 0)
        ret_val = EXIT_FAILURE;

    /* TBD: need to set to H5D_SPARSE_CHUNK */
    if (H5Pset_layout(dcpl, H5D_STRUCT_CHUNK) < 0)
        ret_val = EXIT_FAILURE;

    if (H5Pset_struct_chunk(dcpl, 1, chunk_dim, H5D_SPARSE_CHUNK) < 0)
        ret_val = EXIT_FAILURE;

    if ((did = H5Dcreate2(fid, SPARSE_DSET, H5T_NATIVE_INT, sid, H5P_DEFAULT, dcpl, H5P_DEFAULT)) < 0)
        ret_val = EXIT_FAILURE;

    /* Write sparse data to the dataset */
    // memset(wbuf, 0, sizeof(wbuf));

    /* Add values to the first chunk (for hyperslab selection) */
    fdata[0] = 1;
    fdata[1] = 5;
    fdata[2] = 7;
    
    /*Create the memory space for the first chunk (hyperslab selection)*/
    fmemspace_id = H5Screate_simple(1, fdims, NULL);

    /*Get the dataspace and select the first chunk from the file dataspace*/
    fdataspace_id = H5Dget_space(did);
    H5Sselect_hyperslab(fdataspace_id, H5S_SELECT_SET, foffset, stride, fcount, fblock);

    /*Write the data to the file (via hyperslab selection)*/
    if (H5Dwrite(did, H5T_NATIVE_INT, fmemspace_id, fdataspace_id, H5P_DEFAULT, fdata) < 0)
        ret_val = EXIT_FAILURE;

    printf("The contents of the first subset (chunk), {%d, %d, %d}, written to the file...\n", fdata[0], fdata[1], fdata[2]);

    /*Read from the file (via hyperslab selection)*/
    if (H5Dread(did, H5T_NATIVE_INT, fmemspace_id, fdataspace_id, H5P_DEFAULT, rbuf) < 0 )
        ret_val = EXIT_FAILURE;

    printf("The contents of the first subset (chunk), {%d, %d, %d}, read from the file...\n", rbuf[0], rbuf[1], rbuf[2]);
    
    /* Closing */
    if (H5Sclose(sid) < 0)
        ret_val = EXIT_FAILURE;

    if (H5Pclose(dcpl) < 0)
        ret_val = EXIT_FAILURE;

    if (H5Dclose(did) < 0)
        ret_val = EXIT_FAILURE;

    if (H5Fclose(fid) < 0)
        ret_val = EXIT_FAILURE;

    // PASSED();
    // return SUCCEED;

error:
    H5E_BEGIN_TRY
    {
        H5Sclose(sid);
        H5Sclose(sid1);
        H5Sclose(sid2);
        H5Pclose(dcpl);
        H5Dclose(did);
        H5Fclose(fid);
    }
    H5E_END_TRY

    return FAIL;
}
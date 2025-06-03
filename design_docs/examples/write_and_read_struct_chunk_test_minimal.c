#include "hdf5.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FILE "sparse_write_and_read_chunk_test_minimal.h5"

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
    int rbuf[5];

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
    memset(wbuf, 0, sizeof(wbuf));

    /*ADD IN THE HYPERSLAB SELECTION*/

    /* Initialize and write sparse data to the dataset */
    wbuf[1] = 1;
    wbuf[2] = 3;
    wbuf[3] = 5;
    // wbuf[6] = 6;
    // wbuf[7] = 7;

    /*Write the data to the file */
    if (H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, wbuf) < 0)
        ret_val = EXIT_FAILURE;

    if (H5Dread(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf) < 0 )
        ret_val = EXIT_FAILURE;

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
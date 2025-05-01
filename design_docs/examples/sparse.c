/******************************************************************

  This example shows how to read and write sparse data to a dataset.  
  The data is defined only on a diagonal of the 2D array.

 *****************************************************************/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

#define FILE    "h5_sparse.h5"
#define DATASET "SPARSE_DATA"
#define DIM0    100
#define DIM1    100 
#define CHUNK0   10
#define CHUNK1   10


int
main(void)
{
    hid_t   file   = H5I_INVALID_HID;
    hid_t   space  = H5I_INVALID_HID;
    hid_t   fspace = H5I_INVALID_HID;
    hid_t   dset   = H5I_INVALID_HID;
    hid_t   dcpl   = H5I_INVALID_HID;
    herr_t  status;
    hsize_t dims[2]  = {DIM0, DIM1};
    hsize_t chunk[2] = (CHUNK0, CHUNK1};
     
    int     wdata[DIM0];       /* Write buffer */
    int     rdata[DIM0][DIM1]; /* Read buffer */

    hsize_t      i, j;
    hsize_t      start[2];
    hsize_t      stride[2];
    hsize_t      count[2];
    hsize_t      block[2];

    int num_sections;
    H5_sections_type_t *section_types=NULL;


    /*
     * Initialize data.
     */
    for (i = 0; i < DIM0; i++)
            wdata[i] = i+1 ;

    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple(2, dims, NULL);

    /* 
     * Create a hyperslab selection that corresponds to the diagonal.
     */
     start[0]  = 0;
     start[1]  = 0;
     block[0]  = 1;
     block[1]  = 1;
     stride[0] = 1;
     stride[1] = 1;
     status = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, stride, count, block);
     for (i = 1; i < DIM0; i++) {
         start[0] = i;
         start[1] = i; 
         status = H5Sselect_hyperslab(space, H5S_SELECT_OR, start, stride, count, block);
     }

    /*
     * Create the dataset for storing sparse data.
     */
    dcpl   = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_struct_chunk(dcpl, 2, chunk, H5D_SPARSE_DATA, NULL, NULL); 
    dset   = H5Dcreate(file, DATASET, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write defined data to the dataset.
     */
    status = H5Dwrite(dset, H5T_NATIVE_INT, space, fspace, H5P_DEFAULT, wdata[0]);

    /*
     * Close and release resources.
     */
    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);

    /*
     * Now we begin the read section of this example.
     */

    /*
     * Open file and dataset using the default properties.
     */
    file = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen(file, DATASET, H5P_DEFAULT);

    /*
     * Get information about structured chunk storage.
     */
    dcpl = H5Dget_create_plist(dset);
    if ((H5Pget_layout(dcpl) == H5D_STRUCT_CHUNK) {
        printf ("Layout: H5D_STRUCT_CHUNK\n");
#ifdef 0
        dcpl = H5Pget_struct_chunk_sections(dcpl, &num_sections, NULL);
        if (num_sections > 0) { 
            sections_type = (H5_sections_type_t *)malloc(num_sections * sizeof(H5_sections_type_t));
            dcpl = H5Pget_struct_chunk_sections(dcpl, &num_sections, sections_type); 
            printf("Number of sections in the structired chunk is %d\n", num_sections);
            for (int k; k < num_sections; k++) {
                switch (section_type(k))
                {
                   case H5_SECTION_UNKNOWN:   printf("H5_SECTION_UNKNOWN \n");
                   case H5_SECTION_SELECTION: printf("H5_SECTION_SELECTION\n");
                   case H5_SECTION_FIXED:     printf("H5_SECTION_FIXED\n");             
                   case H5_SECTION_VL:        printf("H5_SECTION_VL\n");             
                   case H5_SECTION_NUM:       printf("Something is not rigth\n");             
                }
             }
#endif
     }
    /*
     * Read the data using the default properties.
     */
    status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata[0]);

    /*
     * Output the data to the screen.
     */
    printf("%s:\n", DATASET);
    for (i = 0; i < DIM0; i++) {
        printf(" [");
        for (j = 0; j < DIM1; j++)
            printf(" %3d", rdata[i][j]);
        printf("]\n");
    }

    /*
     * Close and release resources.
     */
    status = H5Pclose(dcpl);
    status = H5Dclose(dset);
    status = H5Fclose(file);

    return 0;
}

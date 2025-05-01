/******************************************************************

 Example how to read defined elements of the sparse dataset.  

 *****************************************************************/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

#define FILE    "h5_sparse.h5"
#define DATASET "SPARSE_DATA"


int
main(void)
{
    hid_t   file   = H5I_INVALID_HID;
    hid_t   space  = H5I_INVALID_HID;
    hid_t   fspace = H5I_INVALID_HID;
    hid_t   dset   = H5I_INVALID_HID;
    hid_t   dcpl   = H5I_INVALID_HID;
    herr_t  status;
     
    int     *rdata;           /* Read buffer */

    hssize_t num_defined = 0; /*Number of defined elements */

    /*
     * Open file and dataset using the default properties.
     */
    file = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen(file, DATASET, H5P_DEFAULT);

    /*
     * Find and read defined elements. 
     */
    fspace = H5Dget_defined(dset, H5S_ALL, H5P_DEFAULT);

    /*
     * Find number of elements in the selection.
     */
     num_defined = H5Sget_select_npoints(fspace);

     /*
      * Allocate the buffer and read defined elements back.
      */
      rdata = (int *)malloc(num_defined * sizeof(int));
      space = H5Screate_simple(2, dims, NULL);
      status = H5Dread(dset, H5T_NATIVE_INT, space, fspace, H5P_DEFAULT, rdata[0]);

    /*
     * Output the data to the screen.
     */
    printf("Defined elementsi\n");
    for (hssize_t i = 0; i < num_defined; i++) {
            printf(" %3d", rdata[i]);
    }

    /*
     * Close and release resources.
     */
    status = H5Sclose(space);
    status = H5Sclose(fspace);
    status = H5Dclose(dset);
    status = H5Fclose(file);

    return 0;
}

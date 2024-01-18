/*
 * This program generates "n" variable-length (VL) elements of the length between 1 and "m"
 * and stores them in four different HDF5 files:
 *
 * vltype.h5             - contains one 1-dim dataset with generated VL elements using current storage 
 *                         mechanism for variable-length data.
 * vltype_comp.h5        - contains the same dataset but compressed with GZIP deflate level 9.
 * vltype_struct.h5      - contains two 1-dim datasets that emulate storage properties of the structured chunk 
 *                         storage proposed for variable-length data storage: the first dataset contains pairs of 
 *                         offset/length for each element stored in a blob; the second dataset contains the stored
 *                         blob. These two datasets represent two sections of the structured chunk.
 * vltype_struct_comp.h5 - contains the same two datasets but compressed with GZIP deflate level 9.   
 *
 * The program may generate VL vectors with random values or the values that can be successfully compressed
 * based on the value of the flag "d". All datasets use a single chunk to store data.  
 *
 * Use h5dump and h5stat tools to inspect and compare the sizes of stored data and metadata  
 * for the current VL storage approach and for the emulated structured chunk storage.
 * Please remember that for the current VL storage that dataset stores pointers to VL elements; 
 * the elements themselves are stored in the gloabl heap. The h5stat tool reports that space as
 * "Unaccountable space". The h5dump tool -p option will return the size of the dataset with pointers.  
 */

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

#define FILE_NAME1                 		"vltype.h5"
#define FILE_NAME2                 		"vltype_comp.h5"
#define FILE_NAME3                 		"vltype_struct.h5"
#define FILE_NAME4                 		"vltype_struct_comp.h5"
#define VL_DSET_NAME            		"vl_dset"
#define VL_DSET_COMP_NAME	 		"vl_dset_comp"
#define OFFSET_LENGTH_DSET_NAME            	"offset_length_dset"
#define OFFSET_LENGTH_DSET_COMP_NAME 		"offset_length_dset_comp"
#define VL_DATA_DSET_NAME	      		"data"
#define VL_DATA_DSET_COMP_NAME  		"data_comp"
#define NELEMTS     				1000
#define RANK           				1
#define MAX_VL_LEN                      	100

typedef struct {
    long long int   nelemts;
    long long int   max_len;
    int             d;
} handler_t;

handler_t    hand;

/*------------------------------------------------------------
 * Display command line usage
 *------------------------------------------------------------
 */
void
usage(void)
{
    printf("    [-h] [-m --maxLength] [-n --nElements]\n");
    printf("    [-h --help]: this help page\n");
    printf("    [-m --maxLength]: the maximal length of a variable-length element\n");
    printf("    [-n --nElements]: the number of VL type elements in the chunk/dataset\n");
    printf("    [-d --dRandom]: generate random data (default 1) or compressible data (0)\n");
    printf("\n");
}

/*------------------------------------------------------------
 * Parse command line option
 *------------------------------------------------------------
 */
void
parse_command_line(int argc, char *argv[])
{
    int           opt;
    struct option long_options[] = {
                                    {"dimsChunk=", required_argument, NULL, 'c'},
                                    {"help", no_argument, NULL, 'h'},
                                    {"maxLength=", required_argument, NULL, 'm'},
                                    {"nElements=", required_argument, NULL, 'n'},
                                    {"dRandom=", required_argument, NULL, 'd'},
                                    {NULL, 0, NULL, 0}};

    /* Initialize the command line options */
    hand.nelemts = NELEMTS;
    hand.max_len = MAX_VL_LEN;
    hand.d       = 1;
 
    while ((opt = getopt_long(argc, argv, "hm:n:d:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'h':
                fprintf(stdout, "Help page:\n");
                usage();

                exit(0);

                break;
            case 'm':
                /* The maximal length of variable-length element */
                if (optarg) {
                    fprintf(stdout, "maximal length of variable-length element:\t%s\n", optarg);
                    hand.max_len = atoi(optarg);
                }
                else
                    printf("optarg is null\n");
                break;
            case 'n':
                /* The number of VL elements in the dataset */
                if (optarg) {
                    fprintf(stdout, "number of variable-length elements to store:\t%s\n", optarg);
                    hand.nelemts = atoi(optarg);
                }
                else
                    printf("optarg is null\n");
                break;
            case 'd':
                /* The options of data generation */
                if (optarg) {
                    hand.d = atoi(optarg);
                    if (hand.d == 1)
                        fprintf(stdout, "options of data generation:\t\t\trandom values\n");
                    else if (hand.d == 0)
                        fprintf(stdout, "options of data generation:\t\t\tcompressible values\n");
                    else
                        fprintf(stdout, "options of data generation:\t\t\tinvalid option\n");
                }
                else
                    printf("optarg is null\n");
                break;
            case ':':
                printf("option needs a value\n");
                break;
            case '?':
                printf("unknown option: %c\n", optopt);
                break;
        }
    }

    /* optind is for the extra arguments which are not parsed */
    for (; optind < argc; optind++) {
        printf("extra arguments not parsed: %s\n", argv[optind]);
    }

    /* Make sure the command line options are valid */
    if (hand.nelemts < 1 || hand.nelemts > INT_MAX) {
        printf("The number of elements is invalid\n");
        exit(1);
    }

    if (hand.max_len < 1 || hand.max_len > INT_MAX) {
        printf("The maximal length of variable-length element is invalid\n");
        exit(1);
    }

    if (hand.d < 0 || hand.d > 1) {
        printf("Data generation flag can only be 0 (compressible data) or 1 (random)\n");
        exit(1);
    }
}

/*------------------------------------------------------------
 * Create datasets
 *------------------------------------------------------------
 */
int create_dsets(hid_t file, hid_t file_comp, hid_t file_struct, hid_t file_struct_comp)
{
    hid_t   dtype, dcpl, dcpl_compressed;
    hid_t   dset, dset_compressed;
    hid_t   dataspace;
    hvl_t   *vl_data;
    hsize_t dset_dim[1] = {hand.nelemts};
    hsize_t chunk_dim[1] = {hand.nelemts};
    unsigned long long    *the_pairs, *p;
    unsigned long long    offset = 0;
    unsigned long long    total_len = 0;
    char    *all_strings, *ptr;
    int     i, j;

    /* Allocate and initialize variable-length elements */ 
    vl_data = (hvl_t *)malloc(hand.nelemts * sizeof(hvl_t));

    p = the_pairs = (unsigned long long *)malloc(2 * hand.nelemts * sizeof(unsigned long long));

    for(i = 0; i < hand.nelemts; i++) {
        vl_data[i].len = rand() % hand.max_len + 1;
        vl_data[i].p = (char *)malloc(vl_data[i].len * sizeof(char));

        /* Generate random or compressible data */
        for (j = 0; j < vl_data[i].len; j++) {
            if(hand.d)
                ((char *)(vl_data[i].p))[j] = rand() % CHAR_MAX;
            else
            ((char *)(vl_data[i].p))[j] = j;
        }

        /* The pair of offset and length for each variable-length element to be saved in the "structured" dataset */
        *p++ = offset;
        offset += vl_data[i].len;
        *p++ = vl_data[i].len;

        /* The total length of all elements */
        total_len += vl_data[i].len;
    }

    /* Save all VL elements into a buffer for the "structured" dataset */
    ptr = all_strings = (char *)malloc(total_len * sizeof(char));

    for(i = 0; i < hand.nelemts; i++) {
        memcpy(ptr, (char *)(vl_data[i].p), vl_data[i].len);
        ptr += vl_data[i].len;
    }

    dset_dim[0] = hand.nelemts;

    /* Create the dataspace of one chunk size to save the variable-length elements in the current HDF5 way */
    dataspace = H5Screate_simple(RANK, dset_dim, NULL);

    /* Create a VL datatype */
    dtype = H5Tvlen_create(H5T_NATIVE_CHAR);

    dcpl = H5Pcreate(H5P_DATASET_CREATE);

    /* Chunk size must be smaller than 4GB */
    H5Pset_chunk(dcpl, RANK, chunk_dim); 

    /* Create a new dataset without compression to save the variable-length elements in the current HDF5 way */
    dset = H5Dcreate2(file, VL_DSET_NAME, dtype, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    dcpl_compressed = H5Pcopy(dcpl);

    /* Set gzip compression */
    H5Pset_deflate(dcpl_compressed, 9);

    /* Create a new dataset with compression to save the variable-length elements in the current HDF5 way */
    dset_compressed = H5Dcreate2(file_comp, VL_DSET_COMP_NAME, dtype, dataspace, H5P_DEFAULT, dcpl_compressed, H5P_DEFAULT);

    /* Write the data to the dataset */
    H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl_data);

    /* Write the data to the dataset with compression */
    H5Dwrite(dset_compressed, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, vl_data);

    H5Sclose(dataspace);
    H5Dclose(dset);
    H5Dclose(dset_compressed);

    /*--------------------------------------------------------------------------
     * Create and write the dataset for the pairs of offsets and lengths
     * for the proposed "structured" dataset
     *---------------------------------------------------------------------------
     */
    dset_dim[0] = chunk_dim[0] = 2 * hand.nelemts;

    /* Create the dataspace of one chunk size */
    dataspace = H5Screate_simple(RANK, dset_dim, NULL);

    H5Pset_chunk(dcpl, RANK, chunk_dim); 
    H5Pset_chunk(dcpl_compressed, RANK, chunk_dim); 

    /* Create a new dataset without compression */
    dset = H5Dcreate2(file_struct, OFFSET_LENGTH_DSET_NAME, H5T_NATIVE_ULLONG, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    /* Create a new dataset with compression */
    dset_compressed = H5Dcreate2(file_struct_comp, OFFSET_LENGTH_DSET_COMP_NAME, H5T_NATIVE_ULLONG, dataspace, H5P_DEFAULT, dcpl_compressed, H5P_DEFAULT);

    /* Write the data to the dataset */
    H5Dwrite(dset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, the_pairs);

    /* Write the data to the dataset with compression */
    H5Dwrite(dset_compressed, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, the_pairs);

    H5Sclose(dataspace);
    H5Dclose(dset);
    H5Dclose(dset_compressed);

    /*--------------------------------------------------------------------------
     * Create and write the dataset of all VL elements 
     * for the proposed structured dataset
     *---------------------------------------------------------------------------
     */
    dset_dim[0] = chunk_dim[0] = total_len;

    /* Create the dataspace of one chunk size */
    dataspace = H5Screate_simple(RANK, dset_dim, NULL);

    H5Pset_chunk(dcpl, RANK, chunk_dim); 
    H5Pset_chunk(dcpl_compressed, RANK, chunk_dim); 

    /* Create a new dataset without compression */
    dset = H5Dcreate2(file_struct,VL_DATA_DSET_NAME, H5T_NATIVE_CHAR, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    /* Create a new dataset with compression */
    dset_compressed = H5Dcreate2(file_struct_comp, VL_DATA_DSET_COMP_NAME, H5T_NATIVE_CHAR, dataspace, H5P_DEFAULT, dcpl_compressed, H5P_DEFAULT);

    /* Write the data to the dataset */
    H5Dwrite(dset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, all_strings);

    /* Write the data to the dataset with compression */
    H5Dwrite(dset_compressed, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, all_strings);

    H5Sclose(dataspace);
    H5Dclose(dset);
    H5Dclose(dset_compressed);

    /* Free memory buffer */    
    for(i = 0; i < hand.nelemts; i++)
        free(vl_data[i].p);
    free(vl_data);
    free(the_pairs);
    free(all_strings);

    H5Pclose(dcpl_compressed);
    H5Pclose(dcpl);

    return 0;

error:
    return -1;
}

/*------------------------------------------------------------
 * Main function
 *------------------------------------------------------------
 */
int
main(int argc, char **argv)
{
    hid_t   file, file_comp, file_struct, file_struct_comp;
    hsize_t dset_dim[1];
    time_t  t;
    int     n;

    parse_command_line(argc, argv);

    /* Initializing random generator */
    /* srand((unsigned) time(&t)); */

    srand(20);  
    /* Create files */

    file             = H5Fcreate(FILE_NAME1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file_comp        = H5Fcreate(FILE_NAME2, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file_struct      = H5Fcreate(FILE_NAME3, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file_struct_comp = H5Fcreate(FILE_NAME4, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Create datasets */
    create_dsets(file, file_comp, file_struct, file_struct_comp);

    /* Close resources */
    H5Fclose(file);
    H5Fclose(file_comp);
    H5Fclose(file_struct);
    H5Fclose(file_struct_comp);

    return 0;
}

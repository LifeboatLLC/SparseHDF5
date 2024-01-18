/*
 * This program generates a file with a sparse 2-dim data array stored in HDF5 using two options:  
 * as a regular 2-dim HDF5 dataset with a name "sparse" in which undefined elements are represented 
 * by 0 or as two 1-dim HDF5 datasets. In the latter case one dataset with the name "selection"
 * contains the encoded hyperslab selection of the defined elements in the data array and another 
 * dataset with the name "data" contains defined elements themselves. These two 1-dim datasets 
 * emulate two corresponding sections of the structured chunk proposed for sparse storage.
 *
 * All datasets and their compressed counterparts(with the "*_comp" names) are stored under the group 
 * "percent_X", where X indicates percentage of the defined values in the data array and varies between
 * 1 and M (M is specified with a command line option -m). All datasets use chunking storage with 
 * a rank and dimension sizes of the chunk being equal to the rank and dimension sizes of the 
 * corresponding dataset. 
 *
 * The chunk sizes for the 2-dim "sparse" dataset are specified with a command line option -c as multiples of 
 * 1024 bytes. The sizes of the 1-dim HDF5 datasets "data" and "selection" are generated and depend on the i
 * values of X and specified type of the hyperslab selection.
 *  
 * There are three types of hyperslab selection defined by the value specified with a command line option -s:
 *
 *  1 - random locations in each row
 *  2 - random placed rectangular in the entire chunk
 *  3 - randomly placed continuous locations in each row
 *
 * The values of the defined elements of the sparse array are generated based on a command line option -d:
 *
 *  1 - default; the program will generate and initialize data with random values between 1 and UCHAR_MAX
 *  0 - the program will initialize data with the sequences 1,2,...,N, where N =< UCHAR_MAX making it
 *      compressible
 *
 * To compile the program, please use h5cc.  The -h option lists all the command line options:
 *
 *   [-h] [-c --dimsChunk] [-m --mPercent] [-s --spaceSelect] [-d --dRandom] [-v --Verbose]
 * 
 * Example: The commands
 *
 *           h5cc sparse.c
 *           ./a.out -c 1x1 -m 3 -s 2 
 *
 * will generate HDF5 files sparse_file.h5 with the following file structure shown by 
 * the h5dump -n command:
 *
 *	HDF5 "sparse_file.h5" {
 *	FILE_CONTENTS {
 *	 group      /
 *	 group      /percent_1
 *	 dataset    /percent_1/data
 *	 dataset    /percent_1/data_comp
 *	 dataset    /percent_1/selection
 *	 dataset    /percent_1/selection_comp
 *	 dataset    /percent_1/sparse
 *	 dataset    /percent_1/sparse_comp
 *	 group      /percent_2
 *	 dataset    /percent_2/data
 *	 dataset    /percent_2/data_comp
 *	 dataset    /percent_2/selection
 *	 dataset    /percent_2/selection_comp
 *	 dataset    /percent_2/sparse
 *	 dataset    /percent_2/sparse_comp
 *	 group      /percent_3
 *	 dataset    /percent_3/data
 *	 dataset    /percent_3/data_comp
 *	 dataset    /percent_3/selection
 *	 dataset    /percent_3/selection_comp
 *	 dataset    /percent_3/sparse
 *	 dataset    /percent_3/sparse_comp
 *	 }
 *	}
 *
 * Datasets "sparse" and its compressed counterpart "sparse_comp will have 
 * dimensions 1024x1024 (1MiB) with data density varying from 1 to 3%.  All defined elements
 * are located in a randomly placed rectangular sub-region (hyperslab) of the dataset.
 * 
 */

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <getopt.h>

#define FILE_NAME                 	"sparse_file"
#define DSET_NAME	            	"sparse"
#define DSET_COMPRESSED_NAME 		"sparse_comp"
#define DATA_DSET_NAME         		"data"
#define DATA_DSET_COMPRESSED_NAME	"data_comp"
#define SELECTION_DSET_NAME         	"selection"
#define SELECTION_DSET_COMPRESSED_NAME	"selection_comp"
#define GROUP_NAME                	"percent_"
#define GROUP_NUM                 	10
#define CHUNK_DIM1     			10
#define CHUNK_DIM2     			100
#define RANK           			2
#define MAX_PERCENT                     20

typedef struct {
    long long int   chunk_dim1;
    long long int   chunk_dim2;
    int             space_select;
    int             max_percent;
    int             d;               /* flag to generate random or compressible data values */
    int             v;               /* prints progress messages */
} handler_t;

typedef struct {
    long long int   sparse;          /* size of sparse dataset */
    long long int   sparse_comp;     /* size of compressed sparse data set */
    long long int   data;            /* size of dataset with raw data */
    long long int   data_comp;       /* size of comporessed dataste with raw data */
    long long int   sel;             /* size of dataset with encoded selection */
    long long int   sel_comp;        /* size of compressed dataste with encoded selection */
} storage_t;

handler_t    hand;
storage_t    st[MAX_PERCENT];
  

/*------------------------------------------------------------
 * Display command line usage
 *------------------------------------------------------------
 */
void
usage(void)
{
    printf("    [-h] [-c --dimsChunk] [-m --mPercent] [-s --spaceSelect] [-d --dRandom] \n");
    printf("    [-h --help]: this help page\n");
    printf("    [-c --dimsChunk]: the 2D dimensions of the chunks in KB. e.g. 10x20 means the chunk size is 10KB X 20KB.\n");
    printf("    [-m --mPercent]: the maximal percentage of data density, e.g., a value of 5 means the data density will be from 1 to 5 percent.\n");
    printf("	    The datasets will be put into the groups named 'percent_X', where 'X' is 1 to 5. \n");
    printf("    [-s --spaceSelect]: the hyperslab selection of the data density.  The default is random points in each row (value 1).\n");
    printf("	    The other option is an rectangular-shaped selection randomly positioned in the chunk (value 2).\n");
    printf("	    The third option is continuous points in each row with random position (value 3)\n");
    printf("    [-d --dRandom]: Use random data values (1) or compressible data values (0) \n");
    printf("    [-v --Verbose]: Print progress messages(1); default no messages displayed (0) \n");
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
                                    {"mPercent=", required_argument, NULL, 'm'},
                                    {"spaceSelect=", required_argument, NULL, 's'},
                                    {"dRandom=", required_argument, NULL, 'd'},
                                    {"Verbose=", required_argument, NULL, 'v'},
                                    {NULL, 0, NULL, 0}};

    /* Initialize the command line options */
    hand.chunk_dim1               = CHUNK_DIM1; /* First default chunk dimension */
    hand.chunk_dim2               = CHUNK_DIM2; /* Second default chunk dimension */
    hand.space_select             = 1;
    hand.max_percent              = GROUP_NUM;
    hand.d                        = 1;
    hand.v                        = 0;

    while ((opt = getopt_long(argc, argv, "c:hm:s:d:v:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'c':
                /* The dimensions of the chunks */
                if (optarg) {
                    char *dims_str, *dim1_str, *dim2_str;
                    dims_str       = strdup(optarg);
                    dim1_str       = strtok(dims_str, "x");
                    dim2_str       = strtok(NULL, "x");
                    hand.chunk_dim1 = atoll(dim1_str);
                    hand.chunk_dim2 = atoll(dim2_str);
                    /* Input is in KB */
                    hand.chunk_dim1 *= 1024;
                    hand.chunk_dim2 *= 1024;
                    fprintf(stdout, "Chunk dimensions:\t\t\t\t\t\%lld x %lld\n", hand.chunk_dim1, hand.chunk_dim2);
                    free(dims_str);
                }
                else
                    printf("optarg is null\n");
                break;
            case 'h':
                fprintf(stdout, "Help page:\n");
                usage();

                exit(0);

                break;
            case 'm':
                /* The maximal number of percentage of data density, e.g. if given 5, there'll be five groups. 
                 * Each group has a dataset of data density from 1% to 5%. */
                if (optarg) {
                    fprintf(stdout, "Maximal percentage of data density:\t\t\t%s\n", optarg);
                    hand.max_percent = atoi(optarg);
                }
                else
                    printf("optarg is null\n");
                break;
            case 's':
                /* The options of data space selection */
                if (optarg) {
                    hand.space_select = atoi(optarg);

                    if (hand.space_select == 1)
                        fprintf(stdout, "Options of data space selection:\t\t\trandomly selected locations in each row\n");
                    else if (hand.space_select == 2)
                        fprintf(stdout, "Options of data space selection:\t\t\trandomly selected rectangular in the whole chunk\n");
                    else if (hand.space_select == 3)
                        fprintf(stdout, "Options of data space selection:\t\t\trandomly selected continuous locations in each row\n");
                    else
                        fprintf(stdout, "Options of data space selection:\t\t\tinvalid option\n");
                }
                else
                    printf("optarg is null\n");
                break;
           case 'd':
                /* The options of data space selection */
                if (optarg) {
                    hand.d = atoi(optarg);
                    if (hand.d == 1)
                        fprintf(stdout, "Options of data generation:\t\t\t\trandom values\n");
                    else if (hand.d == 0)
                        fprintf(stdout, "Options of data generation:\t\t\t\tcompressible values\n");
                    else
                        fprintf(stdout, "Options of data generation:\t\t\t\tinvalid option\n");
                }
                else
                    printf("optarg is null\n");
                break;
           case 'v':
                /* The options of data space selection */
                if (optarg) {
                    hand.v = atoi(optarg);
                    if (hand.v == 1)
                        fprintf(stdout, "Verbose mode: \t\t\t\t\t\ton\n");
                    else if (hand.v == 0)
                        fprintf(stdout, "Verbose mode: \t\t\t\t\t\toff\n");
                    else
                        fprintf(stdout, "Verbose mode:\t\t\t\t\tinvalid option\n");
                }
                else
                    printf("optarg is null\n");
                break;
            case ':':
                printf("Option needs a value\n");
                break;
            case '?':
                printf("Unknown option: %c\n", optopt);
                break;
        }
    }

    /* optind is for the extra arguments which are not parsed */
    for (; optind < argc; optind++) {
        printf("extra arguments not parsed: %s\n", argv[optind]);
    }

    /* Make sure the command line options are valid */
    if (hand.max_percent < 1 || hand.max_percent > 20) {
        printf("The maximal percentage of the data density isn't valid\n");
        exit(1);
    }

    if (hand.space_select < 1 || hand.space_select > 3) {
        printf("The option of hyperslab selection can only be 1, 2, or 3\n");
        exit(1);
    }

    if (hand.d < 0 || hand.d > 1) {
        printf("Data generation flag can only be 0 (compressible data) or 1 (random)\n");
        exit(1);
    }
    
    if (hand.v < 0 || hand.v > 1) {
        printf("Verbose flag can only be 0 or 1 \n");
        exit(1);
    }
}
/*------------------------------------------------------------
 * Print used storage
 *------------------------------------------------------------
 */
void print_results(int index)
{
   int i;
   long long int a;
   long long int b;
   long long int c;
   float         d;
   

   printf("\n");
   printf("Printing percentage, encoded selection size (ES), compressed encoded selection size (CES), and storage ratio (SR) \n");
   printf("\n");
   printf("         %%         ES        CES         SR\n");
   printf("\n");

   for (i=0; i < index; i++) {
       a = st[i].sel;
       b = st[i].sel_comp;
       d = (float) a/b;
       printf ("%10d %10lli %10lli %10.1f \n", i+1, a, b, d);
   }

   printf("\n");
   printf("Printing percentage, sparse storage size (SPS), structured storage size (STS), and storage ratio (SR) \n");
   printf("\n");
   printf("         %%        SPS        STS         SR\n");
   printf("\n");

   for (i=0; i < index; i++) {
       a = st[i].sparse;
       b = st[i].data;
       c = st[i].sel;
       d = (float)a/(b+c);
       printf ("%10d %10lli %10lli %10.1f \n", i+1, a, b+c, d);
   }

   printf("\n");
   printf("Printing percentage, compressed sparse storage size (CSPS), compressed structured storage size (CSTS), and storage ratio (SR)\n");
   printf("\n");
   printf("         %%       CSPS       CSTS         SR\n");
   printf("\n");

   for (i=0; i < index; i++) { 
       a = st[i].sparse_comp;
       b = st[i].data_comp;
       c = st[i].sel_comp;
       d = (float)a/(b+c);
       printf ("%10d %10lli %10lli %10.1f \n", i+1, a, b+c, d);
   }
   printf("\n");
}    
   

/*------------------------------------------------------------
 * Create compressed and uncompressed datasets to store
 * the encoded dataspace
 *------------------------------------------------------------
 */
int create_encoded_dspace(hid_t group, hid_t dataspace, int index)
{
    hid_t dset, dset_compressed;
    hid_t dspace;
    hid_t dcpl;
    hsize_t dim[1];
    size_t nalloc;
    void  *buf;
    hsize_t offset[1]={0};
    hsize_t chunk_bytes=0;
    hsize_t compressed_chunk_bytes=0;

    H5Sencode(dataspace, NULL, &nalloc, H5P_DEFAULT);
    buf = (void *)malloc(nalloc);

    H5Sencode(dataspace, buf, &nalloc, H5P_DEFAULT);

    dim[0] = nalloc; 
    dcpl = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_chunk(dcpl, 1, dim); 
    dspace = H5Screate_simple(1, dim, NULL);

    /* Create a new dataset without compression */
    dset = H5Dcreate2(group, SELECTION_DSET_NAME, H5T_NATIVE_UCHAR, dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    /* Set gzip compression */
    H5Pset_deflate(dcpl, 9);

    /* Create a new dataset with compression */
    dset_compressed = H5Dcreate2(group, SELECTION_DSET_COMPRESSED_NAME, H5T_NATIVE_UCHAR, dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    /* Write the data to the dataset */
    H5Dwrite(dset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    H5Dget_chunk_storage_size(dset, offset, &chunk_bytes);
    st[index].sel = chunk_bytes;
    

    /* Write the data to the dataset with compression */
    H5Dwrite(dset_compressed, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    H5Dget_chunk_storage_size(dset_compressed, offset, &chunk_bytes);
    st[index].sel_comp = chunk_bytes;
 

    H5Dclose(dset);
    H5Dclose(dset_compressed);
    H5Pclose(dcpl);
    H5Sclose(dspace);

    free(buf);

    return 0;

error:
    return -1;
}

/*------------------------------------------------------------
 * Create compressed and uncompressed datasets to store the 
 * defined data as a one-dimensional array
 *------------------------------------------------------------
 */
int create_structured_dsets(hid_t group, uint64_t nelemts, uint8_t *data, int index)
{
    hid_t dcpl, dataspace;
    hid_t dset, dset_compressed;
    hsize_t dim[1] = {nelemts};
    hsize_t chunk_dim[1] = {nelemts};
    hsize_t offset[1]={0};
    hsize_t chunk_bytes=0;
    hsize_t compressed_chunk_bytes=0;

    dcpl = H5Pcreate(H5P_DATASET_CREATE);

    /* Chunk size must be smaller than 4GB.  Make it slightly smaller than 4GB */
    H5Pset_chunk(dcpl, 1, chunk_dim); 
    
    /* Create dataspace  */
    dataspace = H5Screate_simple(1, dim, NULL);

    /* Create a new dataset without compression */
    dset = H5Dcreate2(group, DATA_DSET_NAME, H5T_STD_U8LE, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    /* Set gzip compression */
    H5Pset_deflate(dcpl, 9);

    /* Create a new dataset without compression */
    dset_compressed = H5Dcreate2(group, DATA_DSET_COMPRESSED_NAME, H5T_STD_U8LE, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);

    /* Write the data to the dataset  and calculate storage*/
    H5Dwrite(dset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dget_chunk_storage_size(dset, offset, &chunk_bytes);
    st[index].data = chunk_bytes;

    /* Write the data to the compressed dataset and calculate storage */
    H5Dwrite(dset_compressed, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dget_chunk_storage_size(dset_compressed, offset, &chunk_bytes);
    st[index].data_comp = chunk_bytes;

    H5Dclose(dset);
    H5Dclose(dset_compressed);
    H5Pclose(dcpl);

    return 0;

error:
    return -1;
}

/*------------------------------------------------------------
 * Create sparse datasets
 *------------------------------------------------------------
 */
int create_hdf5_dsets(hid_t group, hid_t dcpl, hid_t dataspace, uint64_t nelemts, uint8_t *data, int st_index)
{
    char    dset_name[32];
    hid_t   hdf5_dset, hdf5_dset_compressed;   
    hid_t   mem_space;
    hid_t   dcpl_compressed;
    int     d = hand.d;
    hsize_t mem_dim[1];
    hsize_t chunk_offset[2]={0, 0};
    hsize_t chunk_bytes=0;
    hsize_t compressed_chunk_bytes=0;
    herr_t  status;

    /* Create a new dataset without compression */
    hdf5_dset = H5Dcreate2(group, DSET_NAME, H5T_STD_U8LE, dataspace, H5P_DEFAULT, dcpl, H5P_DEFAULT);


    /* Set gzip compression */
    dcpl_compressed = H5Pcopy(dcpl);
    H5Pset_deflate(dcpl_compressed, 9);

    /* Create a new dataset with compression */
    hdf5_dset_compressed = H5Dcreate2(group, DSET_COMPRESSED_NAME, H5T_STD_U8LE, dataspace, H5P_DEFAULT, dcpl_compressed, H5P_DEFAULT);

    mem_dim[0] = nelemts; 
    mem_space = H5Screate_simple(1, mem_dim, NULL);

    /* Write the data to the dataset and calculate storage */
    status = H5Dwrite(hdf5_dset, H5T_NATIVE_UCHAR, mem_space, dataspace, H5P_DEFAULT, data);
    H5Dget_chunk_storage_size(hdf5_dset, chunk_offset, &chunk_bytes);
    st[st_index].sparse = chunk_bytes; 

    /* Write the data to the compressed dataset and calculate storage */
    status = H5Dwrite(hdf5_dset_compressed, H5T_NATIVE_UCHAR, mem_space, dataspace, H5P_DEFAULT, data);
    H5Dget_chunk_storage_size(hdf5_dset_compressed, chunk_offset, &chunk_bytes);
    st[st_index].sparse_comp = chunk_bytes; 

    H5Sclose(mem_space);
    H5Dclose(hdf5_dset);
    H5Dclose(hdf5_dset_compressed);
    H5Pclose(dcpl_compressed);

    return 0;
}

/*------------------------------------------------------------
 * create_hyperslab (dataspace, nelemts);
 *------------------------------------------------------------
 */
uint64_t create_hyperslab(int select_percent, hid_t *dataspace)
{
    int     d = hand.d;
    uint64_t num_selections = hand.chunk_dim2 * select_percent / 100;
    uint64_t sections = 100 / select_percent;
    uint64_t nelemts = 0;
    hsize_t offset[RANK];
    hsize_t block[RANK] = {1, 1};
    int     i, j, k, n;
    int     st_index = (int)select_percent - 1; /* Index of the current st structure element */

    /* The hyperslab selection is defined in three ways:
     *   1. random points in each row.
     *   2. a rectangular randomly positioned in the chunk.
     *   3. continuous points in each row.
     */
    if (hand.space_select == 1) {
        uint64_t num_selections = hand.chunk_dim2 * select_percent / 100;
        uint64_t sections = 100 / select_percent;

        /* Loop through each row and add random points to the hyperslab selection */
        for (i = 0; i < hand.chunk_dim1; i++) {
            offset[0] = i;

            /* Loop through the number of selection (num_selections) according to the selection percentage.
             * There should be one random point being selected in each section. */
            for (j = 0; j < num_selections; j++) {
                offset[1] = j * sections + rand() % sections;

                if (i == 0 && j == 0)
                    H5Sselect_hyperslab(*dataspace, H5S_SELECT_SET, offset, NULL, block, NULL);
                else
                    H5Sselect_hyperslab(*dataspace, H5S_SELECT_OR, offset, NULL, block, NULL);

                /* Total number of points being selected */
                nelemts++;
            }
        }
    } else if (hand.space_select == 2) {
        /* Limit the upper-left corner of the rectangular within the upper-left quadriple of the chunk */
        offset[0] = rand() % (hand.chunk_dim1 / 2);
        offset[1] = rand() % (hand.chunk_dim2 / 2);

        /* Make the rectangular the same shape as the chunk */
        block[0] = hand.chunk_dim1 * sqrt(select_percent) / 10;
        block[1] = hand.chunk_dim2 * sqrt(select_percent) / 10;

        H5Sselect_hyperslab(*dataspace, H5S_SELECT_SET, offset, NULL, block, NULL);

        nelemts = block[0] * block[1];
    } else if (hand.space_select == 3) {
        /* The number of points is fixed to simplify the computation */
        block[0] = 1;
        block[1] = hand.chunk_dim2 * select_percent / 100;

        /* Loop through each row */
        for (i = 0; i < hand.chunk_dim1; i++) {
            /* The position is random in each row */
            offset[0] = i;
            offset[1] = rand() % (hand.chunk_dim2 - block[1]);

            if (i == 0)
                H5Sselect_hyperslab(*dataspace, H5S_SELECT_SET, offset, NULL, block, NULL);
            else
                H5Sselect_hyperslab(*dataspace, H5S_SELECT_OR, offset, NULL, block, NULL);
        }

        /* Total number of points being selected */
        nelemts = hand.chunk_dim1 * block[1];
    }
    return nelemts; 
}

/*------------------------------------------------------------
 * Main function
 *------------------------------------------------------------
 */
int
main(int argc, char **argv)
{
    char    file_name[32];
    char    group_name[32];
    hid_t   file, group;
    hid_t   dcpl, dataspace;
    hsize_t chunk_dims[2];
    time_t  t;
    int     n;
    uint8_t *data, *p;
    uint64_t nelemts = 0;
    uint64_t i;

    parse_command_line(argc, argv);

    /* Initializing random generator; for now use the same seed for reproducibility of the resulst */
    /*    srand((unsigned) time(&t)); */
    srand(2);

    /* Create a new file using H5F_ACC_TRUNC access */
    sprintf(file_name, "%s.h5", FILE_NAME);
    file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dcpl = H5Pcreate(H5P_DATASET_CREATE);

    chunk_dims[0] = hand.chunk_dim1;
    chunk_dims[1] = hand.chunk_dim2;

    /* Chunk size must be smaller than 4GB */
    H5Pset_chunk(dcpl, RANK, chunk_dims); 

    /* Create the dataspace of one chunk size */
    dataspace = H5Screate_simple(RANK, chunk_dims, NULL);

    if (hand.v) printf("Generating file\n");

    for (n = 0; n < hand.max_percent; n++) {

        /* Create some groups in the file. The number is equal to the maximal percentage of the data density */
        sprintf(group_name, "%s%d", GROUP_NAME, n + 1);
        group = H5Gcreate(file, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        /* Generate hyperslab selection and sparse data to store */
        nelemts = create_hyperslab ((n+1), &dataspace);

        /* Generate data */
         p = data = (uint8_t *)malloc(nelemts);

        /* Generate random or compressible values for the defined data */
        for (i = 0; i < nelemts; i++) {
            if (hand.d)
                *p++ = rand() % UCHAR_MAX + 1;
            else
                *p++ = (i+1) % UCHAR_MAX;
        }

        /* Create datasets in the group */
        create_hdf5_dsets(group, dcpl, dataspace, nelemts, data, n);

        /* Create datasets with encoded selection */
        create_encoded_dspace(group, dataspace, n);

        /* Create datasets with defined values */
        create_structured_dsets(group, nelemts, data, n);

        /* Reset hyperslab selection and free data buffer before going to the next iteration*/
        H5Sselect_none(dataspace);
        free(data);

        H5Gclose(group);
        if (hand.v) printf("Closed group %d\n", n+1);
    }

    if (hand.v) printf("Done! \n");

    /* Close resources */
    H5Sclose(dataspace);
    H5Pclose(dcpl);
    H5Fclose(file);

    /* Print results */
    print_results(hand.max_percent);    

    return 0;
}

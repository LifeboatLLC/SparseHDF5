/*
 * mm2h5.c
 *
 * This program reads triplets of unsigned long long integers from a text file
 * and writes them into a two-dimensional chunked and compressed HDF5 dataset.
 *
 * Usage:
 *     ./mm2h5 <file.txt> <XC> <YC> <GZIP_LEVEL> [-v]
 *
 * Behavior:
 *   - The first valid triplet in <file.txt> provides the dataset dimensions (X, Y).
 *   - Each subsequent group of triplets (sharing the same second number) is written
 *     as a hyperslab selection using (a[i], b[i]) as 0-based coordinates and
 *     c[i] as the data value.
 *   - The output HDF5 file is named <file_gzipLEVEL.h5>.
 *   - The dataset is chunked with dimensions (XC, YC) and compressed using GZIP.
 *   - After writing, the program verifies dataset integrity:
 *       • If total points ≤ 10 000 → verifies **all** points
 *       • Otherwise → verifies 10 random points
 *     and prints a verification summary.
 *
 * Options:
 *     -v or --verbose : enable detailed progress output for group writes and verification
 *
 * Example:
 *     ./mm2h5 data.txt 10 10 6 -v
 *     → creates data_gzip6.h5, writes with GZIP=6, prints progress reports
 *
 * Dependencies:
 *     Requires HDF5 C library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "hdf5.h"

#define MAX_LINE_LEN 256
#define VERIFY_POINTS 10
#define VERIFY_ALL_LIMIT 10000

// Reads the next group of triplets where all share the same 'y' value
size_t read_next_group(FILE *fp,
                       unsigned long long **a_out,
                       unsigned long long **b_out,
                       unsigned long long **c_out,
                       unsigned long long *y_value)
{
    char line[MAX_LINE_LEN];
    unsigned long long x, y, z;
    int found = 0;

    // Skip comments and blanks
    while (fgets(line, sizeof(line), fp)) {
        char *ptr = line;
        while (isspace((unsigned char)*ptr)) ptr++;
        if (*ptr == '%' || *ptr == '\0' || *ptr == '\n')
            continue;
        if (sscanf(ptr, "%llu %llu %llu", &x, &y, &z) == 3) {
            found = 1;
            break;
        }
    }
    if (!found)
        return 0;

    size_t capacity = 16, count = 0;
    unsigned long long *a = malloc(capacity * sizeof(unsigned long long));
    unsigned long long *b = malloc(capacity * sizeof(unsigned long long));
    unsigned long long *c = malloc(capacity * sizeof(unsigned long long));
    if (!a || !b || !c) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    unsigned long long current_y = y;
    *y_value = y;

    a[count] = x;
    b[count] = y;
    c[count] = z;
    count++;

    // Continue reading until y changes
    while (fgets(line, sizeof(line), fp)) {
        char *ptr = line;
        while (isspace((unsigned char)*ptr)) ptr++;
        if (*ptr == '%' || *ptr == '\0' || *ptr == '\n')
            continue;
        unsigned long long x2, y2, z2;
        if (sscanf(ptr, "%llu %llu %llu", &x2, &y2, &z2) != 3)
            continue;
        if (y2 != current_y) {
            fseek(fp, -((long)strlen(line)), SEEK_CUR);
            break;
        }
        if (count == capacity) {
            capacity *= 2;
            a = realloc(a, capacity * sizeof(unsigned long long));
            b = realloc(b, capacity * sizeof(unsigned long long));
            c = realloc(c, capacity * sizeof(unsigned long long));
            if (!a || !b || !c) {
                fprintf(stderr, "Reallocation failed\n");
                exit(EXIT_FAILURE);
            }
        }
        a[count] = x2;
        b[count] = y2;
        c[count] = z2;
        count++;
    }

    *a_out = a;
    *b_out = b;
    *c_out = c;
    return count;
}

int main(int argc, char *argv[])
{
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <file.txt> <XC> <YC> <GZIP_LEVEL> [-v]\n", argv[0]);
        return EXIT_FAILURE;
    }

    int verbose = 0;
    if (argc > 5 && (strcmp(argv[5], "-v") == 0 || strcmp(argv[5], "--verbose") == 0))
        verbose = 1;

    const char *txt_file = argv[1];
    hsize_t XC = (hsize_t)atoll(argv[2]);
    hsize_t YC = (hsize_t)atoll(argv[3]);
    int gzip_level = atoi(argv[4]);

    if (gzip_level < 0 || gzip_level > 9) {
        fprintf(stderr, "GZIP level must be between 0 and 9.\n");
        return EXIT_FAILURE;
    }

    // Derive HDF5 filename
    char h5_file[256];
    char base_name[128];
    strncpy(base_name, txt_file, sizeof(base_name));
    char *dot = strrchr(base_name, '.');
    if (dot) *dot = '\0';
    snprintf(h5_file, sizeof(h5_file), "%s_gzip%d.h5", base_name, gzip_level);

    FILE *fp = fopen(txt_file, "r");
    if (!fp) {
        perror("Error opening text file");
        return EXIT_FAILURE;
    }

    // Read first valid triplet (contains dataset dimensions X, Y)
    unsigned long long X = 0, Y = 0, dummy;
    char line[MAX_LINE_LEN];
    while (fgets(line, sizeof(line), fp)) {
        char *ptr = line;
        while (isspace((unsigned char)*ptr)) ptr++;
        if (*ptr == '%' || *ptr == '\0' || *ptr == '\n')
            continue;
        if (sscanf(ptr, "%llu %llu %llu", &X, &Y, &dummy) == 3)
            break;
    }

    if (X == 0 || Y == 0) {
        fprintf(stderr, "Invalid or missing dataset dimensions.\n");
        fclose(fp);
        return EXIT_FAILURE;
    }

    printf("Creating dataset of size %llu x %llu with chunk %llu x %llu, GZIP=%d\n",
           X, Y, XC, YC, gzip_level);

    // Create HDF5 file and dataset
    hid_t file_id = H5Fcreate(h5_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Failed to create HDF5 file.\n");
        fclose(fp);
        return EXIT_FAILURE;
    }

    hsize_t dims[2] = {X, Y};  // (rows, cols)
    hsize_t chunk[2] = {XC, YC};
    hid_t space_id = H5Screate_simple(2, dims, NULL);

    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 2, chunk);
    if (gzip_level > 0)
        H5Pset_deflate(dcpl, gzip_level);

    hid_t dset_id = H5Dcreate2(file_id, "data", H5T_NATIVE_ULLONG, space_id,
                               H5P_DEFAULT, dcpl, H5P_DEFAULT);

    H5Pclose(dcpl);
    H5Sclose(space_id);

    // Store all points for verification
    typedef struct { unsigned long long a, b, c; } point_t;
    point_t *points = NULL;
    size_t total_points = 0, cap_points = 0;
    size_t group_count = 0;

    unsigned long long *a = NULL, *b = NULL, *c = NULL;
    unsigned long long y_value;
    size_t n;

    while ((n = read_next_group(fp, &a, &b, &c, &y_value)) > 0) {
        group_count++;
        for (size_t i = 0; i < n; i++) { a[i]--; b[i]--; }

        // Save points for verification
        if (total_points + n > cap_points) {
            cap_points = (cap_points == 0) ? n * 2 : cap_points * 2;
            points = realloc(points, cap_points * sizeof(point_t));
        }
        for (size_t i = 0; i < n; i++) {
            points[total_points].a = a[i];
            points[total_points].b = b[i];
            points[total_points].c = c[i];
            total_points++;
        }

        hid_t dspace = H5Dget_space(dset_id);
        H5Sselect_none(dspace);
        for (size_t i = 0; i < n; i++) {
            hsize_t start[2] = {a[i], b[i]};
            hsize_t count[2] = {1, 1};
            H5Sselect_hyperslab(dspace, H5S_SELECT_OR, start, NULL, count, NULL);
        }

        hid_t mspace = H5Screate_simple(1, (hsize_t[]){n}, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_ULLONG, mspace, dspace, H5P_DEFAULT, c);
        H5Sclose(mspace);
        H5Sclose(dspace);
        free(a); free(b); free(c);

        if (verbose && group_count % 10 == 0)
            printf("[Write] Processed %zu groups, total %zu points so far...\n",
                   group_count, total_points);
    }

    H5Dclose(dset_id);
    H5Fclose(file_id);
    fclose(fp);

    printf("Data successfully written to %s\n", h5_file);

    // --- Verification ---
    printf("\nStarting verification...\n");
    srand((unsigned)time(NULL));

    file_id = H5Fopen(h5_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset_id = H5Dopen(file_id, "data", H5P_DEFAULT);

    size_t check_count = (total_points <= VERIFY_ALL_LIMIT) ? total_points : VERIFY_POINTS;
    int mismatches = 0;

    if (total_points <= VERIFY_ALL_LIMIT) {
        printf("Verifying all %zu points...\n", total_points);
        for (size_t i = 0; i < total_points; i++) {
            hsize_t start[2] = {points[i].a, points[i].b};
            hsize_t count[2] = {1, 1};
            hid_t fspace = H5Dget_space(dset_id);
            H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
            hid_t mspace = H5Screate_simple(1, (hsize_t[]){1}, NULL);
            unsigned long long value = 0;
            H5Dread(dset_id, H5T_NATIVE_ULLONG, mspace, fspace, H5P_DEFAULT, &value);
            if (value != points[i].c)
                mismatches++;
            if (verbose && i % 1000 == 0 && i > 0)
                printf("[Verify] Checked %zu / %zu points...\n", i, total_points);
            H5Sclose(mspace);
            H5Sclose(fspace);
        }
    } else {
        printf("Dataset has %zu points; verifying %d random samples...\n",
               total_points, VERIFY_POINTS);
        for (int i = 0; i < VERIFY_POINTS; i++) {
            size_t idx = rand() % total_points;
            hsize_t start[2] = {points[idx].a, points[idx].b};
            hsize_t count[2] = {1, 1};
            hid_t fspace = H5Dget_space(dset_id);
            H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, count, NULL);
            hid_t mspace = H5Screate_simple(1, (hsize_t[]){1}, NULL);
            unsigned long long value = 0;
            H5Dread(dset_id, H5T_NATIVE_ULLONG, mspace, fspace, H5P_DEFAULT, &value);
            printf("  (%llu,%llu): expected=%llu, read=%llu %s\n",
                   points[idx].a+1, points[idx].b+1,
                   points[idx].c, value,
                   (value == points[idx].c) ? "OK" : "MISMATCH");
            if (value != points[idx].c)
                mismatches++;
            H5Sclose(mspace);
            H5Sclose(fspace);
        }
    }

    H5Dclose(dset_id);
    H5Fclose(file_id);
    free(points);

    double pct = check_count ? 100.0 * (check_count - mismatches) / check_count : 0.0;
    printf("\nVerification summary:\n");
    printf("  Points checked  : %zu\n", check_count);
    printf("  Mismatches      : %d\n", mismatches);
    printf("  Match rate      : %.2f%%\n", pct);

    printf("Verification complete.\n");
    return EXIT_SUCCESS;
}


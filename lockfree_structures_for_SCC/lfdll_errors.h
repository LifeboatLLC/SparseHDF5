#ifndef LFDLL_ERRORS_H
#define LFDLL_ERRORS_H

#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdatomic.h>
#include <stdlib.h>


/* Default to standalone mode unless HDF5 support is explicitly enabled. */
#ifndef USE_HDF5
#define USE_HDF5 0
#endif

#if USE_HDF5
    /* Only include HDF5 headers when HDF5 support is enabled. */
    #include "H5Epublic.h"
    #include "H5Eprivate.h"
    int err_occurred = FALSE;

    /************************************************************************
    * LFDLL_ERROR()
    *
    * Error macro for the production build of the LFDLL. Uses HDF5's existing
    * error suite. Performs call to HGOTO_ERROR defined in H5Eprivate.h.
    *
    *                                                   AZO -- 5/25/25
    ************************************************************************/
    #define LFDLL_ERROR(maj_id, min_id, ret_val, str)                                                       \
    do {                                                                                                    \
            HGOTO_ERROR(maj_id, min_id, ret_val, str);                                                      \
    } while(0)

#else
    /***************************/
    /* Dummy Major Error Codes */
    /***************************/
    #define H5E_ARGS 0x01
    #define H5E_RESOURCE 0x03

    /***************************/
    /* Dummy Minor Error Codes */
    /***************************/
    #define H5E_BADVALUE 0x02
    #define H5E_CANTALLOC 0x04

    #define SUCCEED 0
    #define FAIL (-1)

    /************************************************************************
    * LFDLL_ERROR()
    *
    * Error macro for the standalone build of the LFDLL. Takes dummy major and
    * minor codes, the return value and a string. Only the string is used here,
    * using the same parameters as the production build allows for an easier
    * and more light-weight invocation. Prints the message/string
    * to the console and asserts() or exits().
    *
    *                                                   AZO -- 5/25/25
    ************************************************************************/
    #define LFDLL_ERROR(maj_id, min_id, ret_val, str)                                                       \
        do {                                                                                                \
                fprintf(stderr,                                                                             \
                    "[LFDLL ERROR] %s:%d: %s \n", __FILE__, __LINE__, str);                                 \
                exit(1);                                                                                    \
                                                                                                            \
        } while(0)


#endif


#endif /* LFDLL_ERRORS_H */

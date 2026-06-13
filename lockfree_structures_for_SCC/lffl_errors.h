#ifndef LFFL_ERRORS_H
#define LFFL_ERRORS_H

#include "lfdll_errors.h"

/************************************************************************
 * LFFL_ERROR()
 *
 * Thin wrapper over LFDLL_ERROR() to keep invocation style consistent.
 ************************************************************************/

#define LFFL_ERROR(maj_id, min_id, ret_val, str)                                                            \
    do {                                                                                                    \
        LFDLL_ERROR((maj_id), (min_id), (ret_val), (str));                                                  \
    } while(0)

#endif // LFFL_ERRORS_H

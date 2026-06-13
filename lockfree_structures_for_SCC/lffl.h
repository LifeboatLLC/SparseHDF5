/************************************************************************
 *
 * File: lffl.h
 *
 *     The lock-free free-list used by the lock-free list macros.
 *
 *     This implementation follows the free-list policy described in the
 *     RFC. The free-list is a lock-free head/tail queue of tokens, with
 *     a stable anchor token present for the lifetime of the queue. The
 *     queue is logically empty when fl_shead == fl_stail and the head
 *     token has no successor.
 *
 *     Each protected operation uses the per-entry token scheme:
 *
 *         - On entry, a thread attempts to reuse the token at the head
 *           of the free-list.
 *
 *         - If that head token still has a non-zero reference count, the
 *           thread leaves it alone and allocates a fresh token instead.
 *
 *         - The selected token is then appended to the tail of the
 *           free-list with ref_count = 1, becoming that thread's pinned
 *           token.
 *
 *         - On exit, the thread decrements the reference count of its
 *           pinned token.
 *
 *     Retired user nodes are published by storing the user node pointer
 *     and deleter in a token and appending that token to the tail of the
 *     free-list with ref_count = 0.
 *
 *     Tokens are only reused from the head of the free-list, and the
 *     links are counted-pointers. A token is not reused while any thread
 *     could still be referencing it, and stale ABA cycles on free-list
 *     links are detected by the serial numbers carried in those links.
 *
 *                                                     CMM -- 03/12/26
 *
 ************************************************************************/
#ifndef LFFL_H
#define LFFL_H

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SUCCEED
#define SUCCEED 0
#endif

#ifndef FAIL
#define FAIL (-1)
#endif

#ifndef LFFL_MALLOC
#define LFFL_MALLOC(sz) malloc(sz)
#endif

#ifndef LFFL_FREE
#define LFFL_FREE(p) free(p)
#endif

#ifndef LFFL_ABORT
#define LFFL_ABORT(msg)                                                                                     \
    do {                                                                                                    \
        fprintf(stderr, "LFFL error at %s:%d: %s\n", __FILE__, __LINE__, (msg));                            \
        abort();                                                                                            \
    } while( 0 )
#endif

#ifndef LFFL_CHECK
#define LFFL_CHECK(cond, msg)                                                                               \
    do {                                                                                                    \
        if(!(cond))                                                                                         \
            LFFL_ABORT(msg);                                                                                \
    } while( 0 )
#endif

#define LFFL_VALID         0x1FF1U
#define LFFL_TOKEN_IN_USE  0x1492U
#define LFFL_TOKEN_ON_FL   0xBEEFU
#define LFFL_TOKEN_INVALID 0xDEADU

typedef void (* lffl_deleter_t)(void * node);

/************************************************************************
 *
 * Counted free-list pointer.
 *
 *     A counted-pointer stores a pointer to a free-list token together
 *     with a serial number. CAS operations compare both fields so that
 *     stale A->B->A link cycles on the free-list do not go undetected.
 *
 ************************************************************************/
typedef struct lffl_sptr_t {
    struct lffl_token_t * ptr;
    unsigned long long    sn;
} lffl_sptr_t;

/************************************************************************
 *
 * Free-list token.
 *
 *     Tokens are the nodes of the lock-free free-list queue. A token may
 *     serve either as a per-entry pinned token or as a token that wraps a
 *     retired user node awaiting safe reclamation.
 *
 *     ref_count is used by the ENTER/EXIT policy described in the RFC.
 *     snext is the counted next link used by the lock-free free-list
 *     queue. node and deleter are NULL when the token is being used only
 *     as an ENTER/EXIT token.
 *
 ************************************************************************/
typedef struct lffl_token_t {
    _Atomic unsigned int  tag;
    _Atomic unsigned int  ref_count;
    _Atomic lffl_sptr_t   snext;

    void                * node;
    lffl_deleter_t        deleter;
} lffl_token_t;

/************************************************************************
 *
 * Statistics block.
 *
 *     These counters are debugging aids only.
 *
 ************************************************************************/
typedef struct lffl_stats_t {
    _Atomic long long enter_calls;
    _Atomic long long exit_calls;
    _Atomic long long retire_calls;

    _Atomic long long num_tokens_allocated;
    _Atomic long long num_tokens_freed;
    _Atomic long long num_tokens_added_to_fl;
    _Atomic long long num_tokens_drawn_from_fl;

    _Atomic long long num_fl_head_update_cols;
    _Atomic long long num_fl_tail_update_cols;
    _Atomic long long num_fl_append_cols;

    _Atomic long long num_fl_req_denied_due_to_empty;
    _Atomic long long num_fl_req_denied_due_to_ref_count;

    _Atomic long long num_token_ref_cnt_incs;
    _Atomic long long num_token_ref_cnt_decs;
} lffl_stats_t;

/************************************************************************
 *
 * LFFL control structure.
 *
 *     The free-list is a lock-free head/tail queue of tokens. fl_shead
 *     and fl_stail are counted-pointers to the queue head and tail. The
 *     queue is never empty; initialization installs a stable anchor token.
 *
 ************************************************************************/
typedef struct lffl_t {
    unsigned int tag;

    _Atomic lffl_sptr_t fl_shead;
    _Atomic lffl_sptr_t fl_stail;
    _Atomic long long   fl_len;

    lffl_stats_t        stats;
} lffl_t;

/************************************************************************
 *
 * LFFL_INIT()
 *
 *     Initialize the supplied LFFL instance. One anchor token is created
 *     and installed as both head and tail.
 *
 ************************************************************************/
#define LFFL_INIT(fl_ptr)                                                                                   \
    do {                                                                                                    \
        lffl_token_t * lffl_init__anchor_ptr = NULL;                                                        \
        lffl_sptr_t    lffl_init__null_sptr  = {NULL, 0ULL};                                                \
        lffl_sptr_t    lffl_init__anchor_sptr;                                                              \
                                                                                                            \
        LFFL_CHECK((fl_ptr) != NULL, "fl_ptr == NULL");                                                     \
                                                                                                            \
        (fl_ptr)->tag = LFFL_VALID;                                                                         \
                                                                                                            \
        atomic_init(&((fl_ptr)->fl_shead), lffl_init__null_sptr);                                           \
        atomic_init(&((fl_ptr)->fl_stail), lffl_init__null_sptr);                                           \
        atomic_init(&((fl_ptr)->fl_len), 1LL);                                                              \
                                                                                                            \
        atomic_init(&((fl_ptr)->stats.enter_calls), 0LL);                                                   \
        atomic_init(&((fl_ptr)->stats.exit_calls), 0LL);                                                    \
        atomic_init(&((fl_ptr)->stats.retire_calls), 0LL);                                                  \
        atomic_init(&((fl_ptr)->stats.num_tokens_allocated), 0LL);                                          \
        atomic_init(&((fl_ptr)->stats.num_tokens_freed), 0LL);                                              \
        atomic_init(&((fl_ptr)->stats.num_tokens_added_to_fl), 0LL);                                        \
        atomic_init(&((fl_ptr)->stats.num_tokens_drawn_from_fl), 0LL);                                      \
        atomic_init(&((fl_ptr)->stats.num_fl_head_update_cols), 0LL);                                       \
        atomic_init(&((fl_ptr)->stats.num_fl_tail_update_cols), 0LL);                                       \
        atomic_init(&((fl_ptr)->stats.num_fl_append_cols), 0LL);                                            \
        atomic_init(&((fl_ptr)->stats.num_fl_req_denied_due_to_empty), 0LL);                                \
        atomic_init(&((fl_ptr)->stats.num_fl_req_denied_due_to_ref_count), 0LL);                            \
        atomic_init(&((fl_ptr)->stats.num_token_ref_cnt_incs), 0LL);                                        \
        atomic_init(&((fl_ptr)->stats.num_token_ref_cnt_decs), 0LL);                                        \
                                                                                                            \
        lffl_init__anchor_ptr = (lffl_token_t *)LFFL_MALLOC(sizeof(*lffl_init__anchor_ptr));                \
        LFFL_CHECK(lffl_init__anchor_ptr != NULL, "malloc failed");                                         \
                                                                                                            \
        atomic_init(&((lffl_init__anchor_ptr)->tag), LFFL_TOKEN_ON_FL);                                     \
        atomic_init(&((lffl_init__anchor_ptr)->ref_count), 0U);                                             \
        atomic_init(&((lffl_init__anchor_ptr)->snext), lffl_init__null_sptr);                               \
        lffl_init__anchor_ptr->node    = NULL;                                                              \
        lffl_init__anchor_ptr->deleter = NULL;                                                              \
                                                                                                            \
        lffl_init__anchor_sptr.ptr = lffl_init__anchor_ptr;                                                 \
        lffl_init__anchor_sptr.sn  = 1ULL;                                                                  \
                                                                                                            \
        atomic_store(&((fl_ptr)->fl_shead), lffl_init__anchor_sptr);                                        \
        atomic_store(&((fl_ptr)->fl_stail), lffl_init__anchor_sptr);                                        \
                                                                                                            \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
        assert(NULL != lffl_init__anchor_ptr);                                                              \
        assert(LFFL_TOKEN_ON_FL == atomic_load(&((lffl_init__anchor_ptr)->tag)));                           \
        assert(0U == atomic_load(&((lffl_init__anchor_ptr)->ref_count)));                                   \
        assert(NULL == atomic_load(&((lffl_init__anchor_ptr)->snext)).ptr);                                 \
        assert(0ULL == atomic_load(&((lffl_init__anchor_ptr)->snext)).sn);                                  \
        assert(NULL == lffl_init__anchor_ptr->node);                                                        \
        assert(NULL == lffl_init__anchor_ptr->deleter);                                                     \
        assert(NULL != atomic_load(&((fl_ptr)->fl_shead)).ptr);                                             \
        assert(NULL != atomic_load(&((fl_ptr)->fl_stail)).ptr);                                             \
        assert(atomic_load(&((fl_ptr)->fl_shead)).ptr == lffl_init__anchor_ptr);                            \
        assert(atomic_load(&((fl_ptr)->fl_stail)).ptr == lffl_init__anchor_ptr);                            \
        assert(atomic_load(&((fl_ptr)->fl_shead)).sn == 1ULL);                                              \
        assert(atomic_load(&((fl_ptr)->fl_stail)).sn == 1ULL);                                              \
        assert(atomic_load(&((fl_ptr)->fl_shead)).ptr == atomic_load(&((fl_ptr)->fl_stail)).ptr);           \
        assert(1LL == atomic_load(&((fl_ptr)->fl_len)));                                                    \
        assert(0LL == atomic_load(&((fl_ptr)->stats.enter_calls)));                                         \
        assert(0LL == atomic_load(&((fl_ptr)->stats.exit_calls)));                                          \
        assert(0LL == atomic_load(&((fl_ptr)->stats.retire_calls)));                                        \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_tokens_allocated)));                                \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_tokens_freed)));                                    \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_tokens_added_to_fl)));                              \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_tokens_drawn_from_fl)));                            \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_fl_head_update_cols)));                             \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_fl_tail_update_cols)));                             \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_fl_append_cols)));                                  \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_fl_req_denied_due_to_empty)));                      \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_fl_req_denied_due_to_ref_count)));                  \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_token_ref_cnt_incs)));                              \
        assert(0LL == atomic_load(&((fl_ptr)->stats.num_token_ref_cnt_decs)));                              \
    } while( 0 ) /* LFFL_INIT() */

/************************************************************************
 *
 * LFFL__RESET_TOKEN()
 *
 *     Reinitialize a token for reuse.
 *
 ************************************************************************/
#define LFFL__RESET_TOKEN(tok_ptr)                                                                          \
    do {                                                                                                    \
        lffl_sptr_t lffl_reset__null_sptr = {NULL, 0ULL};                                                   \
                                                                                                            \
        LFFL_CHECK((tok_ptr) != NULL, "tok_ptr == NULL");                                                   \
                                                                                                            \
        assert(LFFL_TOKEN_IN_USE == atomic_load(&((tok_ptr)->tag)) ||                                       \
               LFFL_TOKEN_ON_FL  == atomic_load(&((tok_ptr)->tag)));                                        \
                                                                                                            \
        atomic_store(&((tok_ptr)->ref_count), 0U);                                                          \
        atomic_store(&((tok_ptr)->snext), lffl_reset__null_sptr);                                           \
        (tok_ptr)->node    = NULL;                                                                          \
        (tok_ptr)->deleter = NULL;                                                                          \
                                                                                                            \
        assert(0U == atomic_load(&((tok_ptr)->ref_count)));                                                 \
        assert(NULL == atomic_load(&((tok_ptr)->snext)).ptr);                                               \
        assert(0ULL == atomic_load(&((tok_ptr)->snext)).sn);                                                \
        assert(NULL == (tok_ptr)->node);                                                                    \
        assert(NULL == (tok_ptr)->deleter);                                                                 \
    } while( 0 ) /* LFFL__RESET_TOKEN() */

/************************************************************************
 *
 * LFFL__CREATE_TOKEN()
 *
 *     Attempt to reuse the token at the head of the free-list. The head
 *     token may only be reused if its reference count is zero. Otherwise,
 *     leave the free-list unchanged and allocate a fresh token.
 *
 *     If the reused token wraps a retired user node, run its deleter
 *     before reinitializing the token for reuse.
 *
 ************************************************************************/
#define LFFL__CREATE_TOKEN(fl_ptr, out_tok_ptr)                                                             \
    do {                                                                                                    \
        bool        lffl_ct__search_done = false;                                                           \
        lffl_sptr_t lffl_ct__shead;                                                                         \
        lffl_sptr_t lffl_ct__stail;                                                                         \
        lffl_sptr_t lffl_ct__snext;                                                                         \
        lffl_sptr_t lffl_ct__test_shead;                                                                    \
        lffl_sptr_t lffl_ct__new_shead;                                                                     \
        lffl_sptr_t lffl_ct__new_stail;                                                                     \
                                                                                                            \
        LFFL_CHECK((fl_ptr) != NULL, "fl_ptr == NULL");                                                     \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
                                                                                                            \
        (out_tok_ptr) = NULL;                                                                               \
                                                                                                            \
        while(!lffl_ct__search_done) {                                                                      \
            lffl_ct__shead      = atomic_load(&((fl_ptr)->fl_shead));                                       \
            lffl_ct__stail      = atomic_load(&((fl_ptr)->fl_stail));                                       \
            lffl_ct__snext      = atomic_load(&((lffl_ct__shead.ptr)->snext));                              \
            lffl_ct__test_shead = atomic_load(&((fl_ptr)->fl_shead));                                       \
                                                                                                            \
            assert(NULL != lffl_ct__shead.ptr);                                                             \
            assert(NULL != lffl_ct__stail.ptr);                                                             \
                                                                                                            \
            if((lffl_ct__test_shead.ptr == lffl_ct__shead.ptr) &&                                           \
               (lffl_ct__test_shead.sn  == lffl_ct__shead.sn)) {                                            \
                                                                                                            \
                if(lffl_ct__shead.ptr == lffl_ct__stail.ptr) {                                              \
                    if(NULL == lffl_ct__snext.ptr) {                                                        \
                        atomic_fetch_add(&((fl_ptr)->stats.num_fl_req_denied_due_to_empty), 1LL);           \
                        lffl_ct__search_done = true;                                                        \
                        break;                                                                              \
                    }                                                                                       \
                                                                                                            \
                    assert(NULL != lffl_ct__snext.ptr);                                                     \
                                                                                                            \
                    lffl_ct__new_stail.ptr = lffl_ct__snext.ptr;                                            \
                    lffl_ct__new_stail.sn  = lffl_ct__stail.sn + 1ULL;                                      \
                                                                                                            \
                    if(!atomic_compare_exchange_strong(&((fl_ptr)->fl_stail),                               \
                                                       &lffl_ct__stail,                                     \
                                                       lffl_ct__new_stail))                                 \
                        atomic_fetch_add(&((fl_ptr)->stats.num_fl_tail_update_cols), 1LL);                  \
                }                                                                                           \
                else {                                                                                      \
                    assert(NULL != lffl_ct__snext.ptr);                                                     \
                                                                                                            \
                    if(atomic_load(&((lffl_ct__shead.ptr)->ref_count)) > 0U) {                              \
                        atomic_fetch_add(&((fl_ptr)->stats.num_fl_req_denied_due_to_ref_count), 1LL);       \
                        lffl_ct__search_done = true;                                                        \
                    }                                                                                       \
                    else {                                                                                  \
                        lffl_ct__new_shead.ptr = lffl_ct__snext.ptr;                                        \
                        lffl_ct__new_shead.sn  = lffl_ct__shead.sn + 1ULL;                                  \
                                                                                                            \
                        if(!atomic_compare_exchange_strong(&((fl_ptr)->fl_shead),                           \
                                                           &lffl_ct__shead,                                 \
                                                           lffl_ct__new_shead)) {                           \
                            atomic_fetch_add(&((fl_ptr)->stats.num_fl_head_update_cols), 1LL);              \
                        }                                                                                   \
                        else {                                                                              \
                            (out_tok_ptr) = lffl_ct__shead.ptr;                                             \
                                                                                                            \
                            assert(NULL != (out_tok_ptr));                                                  \
                            assert(LFFL_TOKEN_ON_FL == atomic_load(&((out_tok_ptr)->tag)));                 \
                            assert(0U == atomic_load(&((out_tok_ptr)->ref_count)));                         \
                                                                                                            \
                            if(NULL != (out_tok_ptr)->node) {                                               \
                                assert(NULL != (out_tok_ptr)->deleter);                                     \
                                ((out_tok_ptr)->deleter)((out_tok_ptr)->node);                              \
                            } else {                                                                        \
                                assert(NULL == (out_tok_ptr)->deleter);                                     \
                            }                                                                               \
                                                                                                            \
                            LFFL__RESET_TOKEN(out_tok_ptr);                                                 \
                            atomic_store(&((out_tok_ptr)->tag), LFFL_TOKEN_IN_USE);                         \
                                                                                                            \
                            assert(LFFL_TOKEN_IN_USE == atomic_load(&((out_tok_ptr)->tag)));                \
                            assert(0U == atomic_load(&((out_tok_ptr)->ref_count)));                         \
                            assert(NULL == (out_tok_ptr)->node);                                            \
                            assert(NULL == (out_tok_ptr)->deleter);                                         \
                                                                                                            \
                            atomic_fetch_sub(&((fl_ptr)->fl_len), 1LL);                                     \
                            atomic_fetch_add(&((fl_ptr)->stats.num_tokens_drawn_from_fl), 1LL);             \
                            assert(0LL <= atomic_load(&((fl_ptr)->fl_len)));                                \
                            lffl_ct__search_done = true;                                                    \
                        }                                                                                   \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
                                                                                                            \
        if(NULL == (out_tok_ptr)) {                                                                         \
            lffl_sptr_t lffl_ct__null_sptr = {NULL, 0ULL};                                                  \
                                                                                                            \
            (out_tok_ptr) = (lffl_token_t *)LFFL_MALLOC(sizeof(*(out_tok_ptr)));                            \
            LFFL_CHECK((out_tok_ptr) != NULL, "malloc failed");                                             \
                                                                                                            \
            atomic_init(&((out_tok_ptr)->tag), LFFL_TOKEN_IN_USE);                                          \
            atomic_init(&((out_tok_ptr)->ref_count), 0U);                                                   \
            atomic_init(&((out_tok_ptr)->snext), lffl_ct__null_sptr);                                       \
            (out_tok_ptr)->node    = NULL;                                                                  \
            (out_tok_ptr)->deleter = NULL;                                                                  \
                                                                                                            \
            atomic_fetch_add(&((fl_ptr)->stats.num_tokens_allocated), 1LL);                                 \
                                                                                                            \
            assert(LFFL_TOKEN_IN_USE == atomic_load(&((out_tok_ptr)->tag)));                                \
            assert(0U == atomic_load(&((out_tok_ptr)->ref_count)));                                         \
            assert(NULL == atomic_load(&((out_tok_ptr)->snext)).ptr);                                       \
            assert(0ULL == atomic_load(&((out_tok_ptr)->snext)).sn);                                        \
            assert(NULL == (out_tok_ptr)->node);                                                            \
            assert(NULL == (out_tok_ptr)->deleter);                                                         \
        }                                                                                                   \
                                                                                                            \
        assert(NULL != (out_tok_ptr));                                                                      \
        assert(LFFL_TOKEN_IN_USE == atomic_load(&((out_tok_ptr)->tag)));                                    \
        assert(0U == atomic_load(&((out_tok_ptr)->ref_count)));                                             \
        assert(NULL == (out_tok_ptr)->node);                                                                \
        assert(NULL == (out_tok_ptr)->deleter);                                                             \
    } while( 0 ) /* LFFL__CREATE_TOKEN() */

/************************************************************************
 *
 * LFFL__DISCARD_TOKEN()
 *
 *     Append the supplied token to the tail of the lock-free free-list
 *     queue.
 *
 ************************************************************************/
#define LFFL__DISCARD_TOKEN(fl_ptr, tok_ptr, expected_ref_count)                                            \
    do {                                                                                                    \
        lffl_sptr_t  lffl_dt__stail;                                                                        \
        lffl_sptr_t  lffl_dt__snext;                                                                        \
        lffl_sptr_t  lffl_dt__test_stail;                                                                   \
        lffl_sptr_t  lffl_dt__new_stail;                                                                    \
        lffl_sptr_t  lffl_dt__new_snext;                                                                    \
        lffl_sptr_t  lffl_dt__null_sptr = {NULL, 0ULL};                                                     \
        unsigned int lffl_dt__in_use_tag = LFFL_TOKEN_IN_USE;                                               \
                                                                                                            \
        LFFL_CHECK((fl_ptr) != NULL, "fl_ptr == NULL");                                                     \
        LFFL_CHECK((tok_ptr) != NULL, "tok_ptr == NULL");                                                   \
                                                                                                            \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
        assert(LFFL_TOKEN_IN_USE == atomic_load(&((tok_ptr)->tag)));                                        \
        assert((expected_ref_count) == atomic_load(&((tok_ptr)->ref_count)));                               \
        assert(((tok_ptr)->node == NULL) == ((tok_ptr)->deleter == NULL));                                  \
                                                                                                            \
        atomic_store(&((tok_ptr)->snext), lffl_dt__null_sptr);                                              \
        assert(NULL == atomic_load(&((tok_ptr)->snext)).ptr);                                               \
                                                                                                            \
        if(!atomic_compare_exchange_strong(&((tok_ptr)->tag),                                               \
                                           &lffl_dt__in_use_tag,                                            \
                                           LFFL_TOKEN_ON_FL))                                               \
            assert(false);                                                                                  \
                                                                                                            \
        assert(LFFL_TOKEN_ON_FL == atomic_load(&((tok_ptr)->tag)));                                         \
                                                                                                            \
        while(true) {                                                                                       \
            lffl_dt__stail      = atomic_load(&((fl_ptr)->fl_stail));                                       \
            lffl_dt__snext      = atomic_load(&((lffl_dt__stail.ptr)->snext));                              \
            lffl_dt__test_stail = atomic_load(&((fl_ptr)->fl_stail));                                       \
                                                                                                            \
            assert(NULL != lffl_dt__stail.ptr);                                                             \
                                                                                                            \
            if((lffl_dt__test_stail.ptr == lffl_dt__stail.ptr) &&                                           \
               (lffl_dt__test_stail.sn  == lffl_dt__stail.sn)) {                                            \
                                                                                                            \
                if(NULL == lffl_dt__snext.ptr) {                                                            \
                    lffl_dt__new_snext.ptr = (tok_ptr);                                                     \
                    lffl_dt__new_snext.sn  = lffl_dt__snext.sn + 1ULL;                                      \
                                                                                                            \
                    if(atomic_compare_exchange_strong(&((lffl_dt__stail.ptr)->snext),                       \
                                                      &lffl_dt__snext,                                      \
                                                      lffl_dt__new_snext)) {                                \
                        lffl_dt__new_stail.ptr = (tok_ptr);                                                 \
                        lffl_dt__new_stail.sn  = lffl_dt__stail.sn + 1ULL;                                  \
                                                                                                            \
                        if(!atomic_compare_exchange_strong(&((fl_ptr)->fl_stail),                           \
                                                           &lffl_dt__stail,                                 \
                                                           lffl_dt__new_stail))                             \
                            atomic_fetch_add(&((fl_ptr)->stats.num_fl_tail_update_cols), 1LL);              \
                                                                                                            \
                        atomic_fetch_add(&((fl_ptr)->fl_len), 1LL);                                         \
                        atomic_fetch_add(&((fl_ptr)->stats.num_tokens_added_to_fl), 1LL);                   \
                                                                                                            \
                        assert(1LL <= atomic_load(&((fl_ptr)->fl_len)));                                    \
                        assert(LFFL_TOKEN_ON_FL == atomic_load(&((tok_ptr)->tag)));                         \
                        assert((expected_ref_count) == atomic_load(&((tok_ptr)->ref_count)));               \
                        break;                                                                              \
                    }                                                                                       \
                                                                                                            \
                    atomic_fetch_add(&((fl_ptr)->stats.num_fl_append_cols), 1LL);                           \
                }                                                                                           \
                else {                                                                                      \
                                                                                                            \
                    lffl_dt__new_stail.ptr = lffl_dt__snext.ptr;                                            \
                    lffl_dt__new_stail.sn  = lffl_dt__stail.sn + 1ULL;                                      \
                                                                                                            \
                    if(!atomic_compare_exchange_strong(&((fl_ptr)->fl_stail),                               \
                                                       &lffl_dt__stail,                                     \
                                                       lffl_dt__new_stail))                                 \
                        atomic_fetch_add(&((fl_ptr)->stats.num_fl_tail_update_cols), 1LL);                  \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
    } while( 0 ) /* LFFL__DISCARD_TOKEN() */

/************************************************************************
 *
 * LFFL_ENTER()
 *
 *     Obtain a token for the current operation. The thread first tries
 *     to reuse the token at the head of the free-list. If that head
 *     token still has a positive reference count, a fresh token is
 *     allocated instead. The selected token is then appended to the
 *     tail of the free-list with ref_count = 1 and returned as the
 *     caller's pinned token.
 *
 ************************************************************************/
#define LFFL_ENTER(fl_ptr, tok_ptr_ptr)                                                                     \
    do {                                                                                                    \
        lffl_token_t * lffl_enter__tok_ptr = NULL;                                                          \
                                                                                                            \
        LFFL_CHECK((fl_ptr) != NULL, "fl_ptr == NULL");                                                     \
        LFFL_CHECK((tok_ptr_ptr) != NULL, "tok_ptr_ptr == NULL");                                           \
                                                                                                            \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
                                                                                                            \
        atomic_fetch_add(&((fl_ptr)->stats.enter_calls), 1LL);                                              \
                                                                                                            \
        LFFL__CREATE_TOKEN((fl_ptr), lffl_enter__tok_ptr);                                                  \
                                                                                                            \
        assert(NULL != lffl_enter__tok_ptr);                                                                \
        assert(0U == atomic_load(&((lffl_enter__tok_ptr)->ref_count)));                                     \
        assert(NULL == lffl_enter__tok_ptr->node);                                                          \
        assert(NULL == lffl_enter__tok_ptr->deleter);                                                       \
        assert(LFFL_TOKEN_IN_USE == atomic_load(&((lffl_enter__tok_ptr)->tag)));                            \
                                                                                                            \
        atomic_store(&((lffl_enter__tok_ptr)->ref_count), 1U);                                              \
        atomic_fetch_add(&((fl_ptr)->stats.num_token_ref_cnt_incs), 1LL);                                   \
                                                                                                            \
        assert(1U == atomic_load(&((lffl_enter__tok_ptr)->ref_count)));                                     \
                                                                                                            \
        LFFL__DISCARD_TOKEN((fl_ptr), lffl_enter__tok_ptr, 1U);                                             \
                                                                                                            \
        assert(LFFL_TOKEN_ON_FL == atomic_load(&((lffl_enter__tok_ptr)->tag)));                             \
        assert(1U == atomic_load(&((lffl_enter__tok_ptr)->ref_count)));                                     \
        assert(NULL == lffl_enter__tok_ptr->node);                                                          \
        assert(NULL == lffl_enter__tok_ptr->deleter);                                                       \
                                                                                                            \
        *(tok_ptr_ptr) = lffl_enter__tok_ptr;                                                               \
        assert(NULL != *(tok_ptr_ptr));                                                                     \
        assert(*(tok_ptr_ptr) == lffl_enter__tok_ptr);                                                      \
    } while( 0 ) /* LFFL_ENTER() */

/************************************************************************
 *
 * LFFL_RETIRE()
 *
 *     Publish a retired user node by storing the node pointer and deleter
 *     in a token and appending that token to the tail of the free-list
 *     with ref_count = 0.
 *
 ************************************************************************/
#define LFFL_RETIRE(fl_ptr, node_ptr, deleter_fn)                                                           \
    do {                                                                                                    \
        lffl_token_t * lffl_retire__tok_ptr = NULL;                                                         \
                                                                                                            \
        LFFL_CHECK((fl_ptr) != NULL, "fl_ptr == NULL");                                                     \
        LFFL_CHECK((node_ptr) != NULL, "node_ptr == NULL");                                                 \
        LFFL_CHECK((deleter_fn) != NULL, "deleter_fn == NULL");                                             \
                                                                                                            \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
                                                                                                            \
        atomic_fetch_add(&((fl_ptr)->stats.retire_calls), 1LL);                                             \
                                                                                                            \
        LFFL__CREATE_TOKEN((fl_ptr), lffl_retire__tok_ptr);                                                 \
                                                                                                            \
        assert(NULL != lffl_retire__tok_ptr);                                                               \
        assert(0U == atomic_load(&((lffl_retire__tok_ptr)->ref_count)));                                    \
        assert(NULL == lffl_retire__tok_ptr->node);                                                         \
        assert(NULL == lffl_retire__tok_ptr->deleter);                                                      \
        assert(LFFL_TOKEN_IN_USE == atomic_load(&((lffl_retire__tok_ptr)->tag)));                           \
                                                                                                            \
        lffl_retire__tok_ptr->node    = (node_ptr);                                                         \
        lffl_retire__tok_ptr->deleter = (deleter_fn);                                                       \
                                                                                                            \
        assert(NULL != lffl_retire__tok_ptr->node);                                                         \
        assert(NULL != lffl_retire__tok_ptr->deleter);                                                      \
        assert(0U == atomic_load(&((lffl_retire__tok_ptr)->ref_count)));                                    \
                                                                                                            \
        LFFL__DISCARD_TOKEN((fl_ptr), lffl_retire__tok_ptr, 0U);                                            \
                                                                                                            \
        assert(LFFL_TOKEN_ON_FL == atomic_load(&((lffl_retire__tok_ptr)->tag)));                            \
        assert(0U == atomic_load(&((lffl_retire__tok_ptr)->ref_count)));                                    \
        assert(NULL != lffl_retire__tok_ptr->node);                                                         \
        assert(NULL != lffl_retire__tok_ptr->deleter);                                                      \
    } while( 0 ) /* LFFL_RETIRE() */

/************************************************************************
 *
 * LFFL_EXIT()
 *
 *     Release the caller's pinned token by decrementing its reference
 *     count.
 *
 ************************************************************************/
#define LFFL_EXIT(fl_ptr, tok_ptr)                                                                          \
    do {                                                                                                    \
        unsigned int lffl_exit__old_ref_count;                                                              \
                                                                                                            \
        LFFL_CHECK((fl_ptr) != NULL, "fl_ptr == NULL");                                                     \
        LFFL_CHECK((tok_ptr) != NULL, "tok_ptr == NULL");                                                   \
                                                                                                            \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
        assert(LFFL_TOKEN_ON_FL == atomic_load(&((tok_ptr)->tag)));                                         \
        assert(atomic_load(&((tok_ptr)->ref_count)) > 0U);                                                  \
        assert(NULL == (tok_ptr)->node);                                                                    \
        assert(NULL == (tok_ptr)->deleter);                                                                 \
                                                                                                            \
        atomic_fetch_add(&((fl_ptr)->stats.exit_calls), 1LL);                                               \
                                                                                                            \
        lffl_exit__old_ref_count = atomic_fetch_sub(&((tok_ptr)->ref_count), 1U);                           \
        assert(lffl_exit__old_ref_count > 0U);                                                              \
                                                                                                            \
        atomic_fetch_add(&((fl_ptr)->stats.num_token_ref_cnt_decs), 1LL);                                   \
    } while( 0 ) /* LFFL_EXIT() */

/************************************************************************
 *
 * LFFL_DUMP_STATS()
 *
 *     Dump debugging counters to file_ptr.
 *
 ************************************************************************/
#define LFFL_DUMP_STATS(fl_ptr, file_ptr)                                                                   \
    do {                                                                                                    \
        FILE * lffl_dump__file_ptr = (file_ptr);                                                            \
                                                                                                            \
        LFFL_CHECK((fl_ptr) != NULL, "fl_ptr == NULL");                                                     \
        LFFL_CHECK(lffl_dump__file_ptr != NULL, "file_ptr == NULL");                                        \
                                                                                                            \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
        assert(0LL <= atomic_load(&((fl_ptr)->fl_len)));                                                    \
                                                                                                            \
        fprintf(lffl_dump__file_ptr,                                                                        \
                "LFFL stats: enter=%lld exit=%lld retire=%lld alloc=%lld freed=%lld on_fl=%lld "            \
                "drawn=%lld append_col=%lld head_col=%lld tail_col=%lld empty_denied=%lld "                 \
                "ref_denied=%lld ref_inc=%lld ref_dec=%lld len=%lld \n",                                    \
                (long long)atomic_load(&((fl_ptr)->stats.enter_calls)),                                     \
                (long long)atomic_load(&((fl_ptr)->stats.exit_calls)),                                      \
                (long long)atomic_load(&((fl_ptr)->stats.retire_calls)),                                    \
                (long long)atomic_load(&((fl_ptr)->stats.num_tokens_allocated)),                            \
                (long long)atomic_load(&((fl_ptr)->stats.num_tokens_freed)),                                \
                (long long)atomic_load(&((fl_ptr)->stats.num_tokens_added_to_fl)),                          \
                (long long)atomic_load(&((fl_ptr)->stats.num_tokens_drawn_from_fl)),                        \
                (long long)atomic_load(&((fl_ptr)->stats.num_fl_append_cols)),                              \
                (long long)atomic_load(&((fl_ptr)->stats.num_fl_head_update_cols)),                         \
                (long long)atomic_load(&((fl_ptr)->stats.num_fl_tail_update_cols)),                         \
                (long long)atomic_load(&((fl_ptr)->stats.num_fl_req_denied_due_to_empty)),                  \
                (long long)atomic_load(&((fl_ptr)->stats.num_fl_req_denied_due_to_ref_count)),              \
                (long long)atomic_load(&((fl_ptr)->stats.num_token_ref_cnt_incs)),                          \
                (long long)atomic_load(&((fl_ptr)->stats.num_token_ref_cnt_decs)),                          \
                (long long)atomic_load(&((fl_ptr)->fl_len)));                                               \
    } while( 0 ) /* LFFL_DUMP_STATS() */

/************************************************************************
 *
 * LFFL_DESTROY()
 *
 *     Destroy the free-list after all concurrent access has stopped.
 *
 *     This walks the free-list queue from head to tail, reclaims any
 *     retired user nodes still wrapped by tokens, frees all tokens,
 *     and invalidates the control structure.
 *
 ************************************************************************/
#define LFFL_DESTROY(fl_ptr)                                                                                \
    do {                                                                                                    \
        lffl_token_t * lffl_destroy__tok_ptr;                                                               \
        lffl_token_t * lffl_destroy__next_ptr;                                                              \
        void * lffl_destroy__fl_vp = (void *)(fl_ptr);                                                      \
                                                                                                            \
        if(!lffl_destroy__fl_vp)                                                                            \
            LFFL_ABORT("fl_ptr == NULL");                                                                   \
                                                                                                            \
        assert(LFFL_VALID == (fl_ptr)->tag);                                                                \
                                                                                                            \
        lffl_destroy__tok_ptr = atomic_load(&((fl_ptr)->fl_shead)).ptr;                                     \
        assert(NULL != lffl_destroy__tok_ptr);                                                              \
                                                                                                            \
        while(NULL != lffl_destroy__tok_ptr) {                                                              \
            lffl_destroy__next_ptr = atomic_load(&((lffl_destroy__tok_ptr)->snext)).ptr;                    \
                                                                                                            \
            assert(LFFL_TOKEN_ON_FL == atomic_load(&((lffl_destroy__tok_ptr)->tag)));                       \
            assert(0U == atomic_load(&((lffl_destroy__tok_ptr)->ref_count)));                               \
                                                                                                            \
            if(NULL != lffl_destroy__tok_ptr->node) {                                                       \
                assert(NULL != lffl_destroy__tok_ptr->deleter);                                             \
                (lffl_destroy__tok_ptr->deleter)(lffl_destroy__tok_ptr->node);                              \
            } else {                                                                                        \
                assert(NULL == lffl_destroy__tok_ptr->deleter);                                             \
            }                                                                                               \
                                                                                                            \
            atomic_store(&((lffl_destroy__tok_ptr)->tag), LFFL_TOKEN_INVALID);                              \
            LFFL_FREE(lffl_destroy__tok_ptr);                                                               \
            atomic_fetch_add(&((fl_ptr)->stats.num_tokens_freed), 1LL);                                     \
                                                                                                            \
            lffl_destroy__tok_ptr = lffl_destroy__next_ptr;                                                 \
        }                                                                                                   \
                                                                                                            \
        (fl_ptr)->tag = 0U;                                                                                 \
        atomic_store(&((fl_ptr)->fl_len), 0LL);                                                             \
                                                                                                            \
        assert(0U == (fl_ptr)->tag);                                                                        \
        assert(0LL == atomic_load(&((fl_ptr)->fl_len)));                                                    \
    } while( 0 ) /* LFFL_DESTROY() */

#endif /* LFFL_H */
/******************************************************************************
 * lfdll.h
 *
 * Lock-Free Doubly Linked List (Sundell-Tsigas)
 ******************************************************************************/

#ifndef LFDLL_H
#define LFDLL_H

#include "lfdll_errors.h"
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef LFDLL_MALLOC
#define LFDLL_MALLOC(sz) malloc(sz)
#endif

#ifndef LFDLL_FREE
#define LFDLL_FREE(p) free(p)
#endif

#include "lffl.h"

/* PAUSE and YIELD were previously used in the same spirit as the back-off */
/* steps in the Sundell-Tsigas algorithm. They are not strictly required   */
/* to avoid livelock, but they may provide benefits under contention.      */
/* As such, they are retained as optional behavior.                        */

#ifndef LFDLL_ENABLE_BACKOFF
#define LFDLL_ENABLE_BACKOFF 0
#endif

#if LFDLL_ENABLE_BACKOFF
#include <sched.h>
#if defined(__i386__) || defined(__x86_64__)
#define LFDLL_PAUSE() do { __asm__ __volatile__("pause" ::: "memory"); } while( 0 )
#else
#define LFDLL_PAUSE() do { } while( 0 )
#endif
#define LFDLL_YIELD() do { sched_yield(); } while( 0 )
#else
#define LFDLL_PAUSE() do { } while( 0 )
#define LFDLL_YIELD() do { } while( 0 )
#endif
#define LFDLL_BACKOFF() do { LFDLL_PAUSE(); LFDLL_YIELD(); } while( 0 )

/************************************************************************
 * Safeguard bounds for PREV corrections.
 *                                           
 * LFDLL_CORRECT_PREV_MAX_SPINS limits repair attempts before restarting 
 * the walk. LFDLL_CORRECT_PREV_MAX_HOPS limits traversal distance during 
 * correction/restart. These tune recovery behavior, they do not limit
 * the size of the DLL.
 ************************************************************************/
#ifndef LFDLL_CORRECT_PREV_MAX_SPINS
#define LFDLL_CORRECT_PREV_MAX_SPINS 8192U
#endif

#ifndef LFDLL_CORRECT_PREV_MAX_HOPS
#define LFDLL_CORRECT_PREV_MAX_HOPS 200000U
#endif

/************************************************************************
 *
 * Removal wait/scan tuning parameters.
 *
 *     In the deleter == NULL test configuration, the remove path performs a
 *     bounded helping/verification pass before reporting success. These
 *     limits bound that extra cleanup work without changing the underlying
 *     algorithm.
 *
 ************************************************************************/
#ifndef LFDLL_REMOVE_WAIT_MAX_ROUNDS
#define LFDLL_REMOVE_WAIT_MAX_ROUNDS 128U
#endif

#ifndef LFDLL_REMOVE_WAIT_SCAN_LIMIT
#define LFDLL_REMOVE_WAIT_SCAN_LIMIT 100000U
#endif


struct lfdll_list_t;
typedef struct lfdll_list_t lfdll_list_t;


/************************************************************************
 *
 * Lock-Free Doubly-Linked List
 *
 *     Macro implementation of the Sundell-Tsigas lock-free doubly
 *     linked list. The list is treated as a singly linked list with
 *     auxiliary predecessor information in each node. Forward links are
 *     authoritative, backward links are repaired by helping, and logical
 *     deletion is represented by storing the deletion mark in the low bit
 *     of the link pointer.
 *
 *     The API is made typeless by using caller-supplied link offsets.
 *     The core logic follows the Sundell-Tsigas algorithms found in
 *     both the thesis and the presentation.
 *
 *                                             CMM -- 3/4/26
 *
 ************************************************************************/



/************************************************************************
 *
 * Tag / mark helpers.
 *
 *     The deletion mark is stored in the low bit of each link pointer.
 *
 ************************************************************************/

#define LFDLL_TAG_MASK        ((uintptr_t)1u)
#define LFDLL_IS_MARKED(p)    (((uintptr_t)(p) & LFDLL_TAG_MASK) != 0)
#define LFDLL_MARK(p)         ((uintptr_t)(p) | LFDLL_TAG_MASK)
#define LFDLL_UNMARK(p)       ((uintptr_t)(p) & ~LFDLL_TAG_MASK)
#define LFDLL_GET_PTR(p)      ((void *)LFDLL_UNMARK(p))

#ifndef LFDLL_PTR
#define LFDLL_PTR(p)          LFDLL_GET_PTR(p)
#endif

#define LFDLL_VALID 0x4c46444c


/************************************************************************
 *
 * Internal sentinel and list data structures.
 *
 *     head_s and tail_s are permanent sentinels. The caller supplies the
 *     offsets of the embedded next and prev link fields in the concrete
 *     node type, allowing the list logic to remain type-agnostic.
 *
 ************************************************************************/

typedef struct lfdll_sentinel_t {
    _Atomic(uintptr_t) next;   /* tagged pointer */
    _Atomic(uintptr_t) prev;   /* tagged pointer */
} lfdll_sentinel_t;

struct lfdll_list_t {
    int tag; /* validity tag */

    /* offsets into user node type */
    size_t off_next;
    size_t off_prev;

    /* head/tail sentinels */
    lfdll_sentinel_t head_s;
    lfdll_sentinel_t tail_s;

    /* free-list */
    lffl_t fl;

    /* stats */
    _Atomic long long ins_ok;
    _Atomic long long ins_fail;
    _Atomic long long rem_ok;
    _Atomic long long rem_fail;
    _Atomic long long correct_prev_calls;
    _Atomic long long correct_prev_iters;
};

/************************************************************************
 *
 * Internal helper macros.
 *
 *     These helpers correspond directly to the primitive operations used by
 *     the Sundell-Tsigas style algorithm: load a tagged link, extract the
 *     underlying pointer, mark a link for logical deletion, find/correct a
 *     predecessor, and complete insertion helping.
 *
 ************************************************************************/
/************************************************************************
 *
 * LFDLL__LINK_NEXT()
 *
 *     Return the address of the atomic next link associated with node_ptr.
 *
 *     node_ptr may reference either a user node or one of the two internal
 *     sentinels. The macro uses the stored next-field offset for user nodes
 *     and the sentinel storage for head_s/tail_s.
 *
 ************************************************************************/
#define LFDLL__LINK_NEXT(lst_ptr, node_ptr)                                                                 \
    (                                                                                                       \
        assert(lst_ptr),                                                                                    \
        ((node_ptr) == (void *)&((lst_ptr)->head_s)) ? &((lst_ptr)->head_s.next) :                          \
        ((node_ptr) == (void *)&((lst_ptr)->tail_s)) ? &((lst_ptr)->tail_s.next) :                          \
        (assert(node_ptr), (_Atomic(uintptr_t) *)((char *)(node_ptr) + (lst_ptr)->off_next))                \
    )

/************************************************************************
 *
 * LFDLL__LINK_PREV()
 *
 *     Return the address of the atomic prev link associated with node_ptr.
 *
 *     As with LFDLL__LINK_NEXT(), this macro handles both user nodes and the
 *     internal sentinels. Backward links are advisory and may temporarily lag
 *     until repaired by helping operations.
 *
 ************************************************************************/
#define LFDLL__LINK_PREV(lst_ptr, node_ptr)                                                                 \
    (                                                                                                       \
        assert(lst_ptr),                                                                                    \
        ((node_ptr) == (void *)&((lst_ptr)->head_s)) ? &((lst_ptr)->head_s.prev) :                          \
        ((node_ptr) == (void *)&((lst_ptr)->tail_s)) ? &((lst_ptr)->tail_s.prev) :                          \
        (assert(node_ptr), (_Atomic(uintptr_t) *)((char *)(node_ptr) + (lst_ptr)->off_prev))                \
    )

/************************************************************************
 *
 * LFDLL__LOAD_LINK()
 *
 *     Load a raw link value and report whether the deletion mark is set.
 *
 *     The macro reads either the next or prev link, returns the raw tagged
 *     uintptr_t value through raw_out, and returns the logical deletion state
 *     through marked_out.
 *
 ************************************************************************/
#define LFDLL__LOAD_LINK(lst_ptr, node_ptr, is_next, marked_out, raw_out)                                   \
    do {                                                                                                    \
        _Atomic(uintptr_t) * lfdll_load_link__link_ptr = (is_next) ?                                        \
            LFDLL__LINK_NEXT((lst_ptr), (node_ptr)) :                                                       \
            LFDLL__LINK_PREV((lst_ptr), (node_ptr));                                                        \
        uintptr_t lfdll_load_link__v = atomic_load(lfdll_load_link__link_ptr);                              \
        *(marked_out) = LFDLL_IS_MARKED(lfdll_load_link__v);                                                \
        *(raw_out) = lfdll_load_link__v;                                                                    \
    } while( 0 )

/************************************************************************
 *
 * LFDLL__LOAD_PTR()
 *
 *     Load a link and return only the underlying unmarked pointer.
 *
 *     This is the common helper used by traversal code when the pointer value
 *     is needed without the logical deletion mark.
 *
 ************************************************************************/
#define LFDLL__LOAD_PTR(lst_ptr, node_ptr, is_next, ptr_out)                                                \
    do {                                                                                                    \
        bool lfdll_load_ptr__marked;                                                                        \
        uintptr_t lfdll_load_ptr__raw;                                                                      \
        LFDLL__LOAD_LINK((lst_ptr), (node_ptr), (is_next), &lfdll_load_ptr__marked, &lfdll_load_ptr__raw);  \
        *(ptr_out) = LFDLL_GET_PTR(lfdll_load_ptr__raw);                                                    \
    } while( 0 )

/************************************************************************
 *
 * LFDLL__CAS_LINK()
 *
 *     Perform a strong CAS operation on a tagged link field.
 *     This helper simply keeps all link CAS operations uniform.
 *
 ************************************************************************/
#define LFDLL__CAS_LINK(link_ptr, expect_val, desired_val)                                                  \
    atomic_compare_exchange_strong((link_ptr), &(uintptr_t){(expect_val)}, (desired_val))


/************************************************************************
 *
 * LFDLL__SET_MARK()
 *
 *     Atomically set the deletion mark on either the next or prev link.
 *
 *     This is the local realization of the SetMark helper from the reference
 *     algorithm. The macro loops until the target link is observed marked or
 *     until the marking CAS succeeds.
 *
 ************************************************************************/
#define LFDLL__SET_MARK(lst_ptr, node_ptr, is_next)                                                         \
    do {                                                                                                    \
        assert(lst_ptr);                                                                                    \
        assert(node_ptr);                                                                                   \
        _Atomic(uintptr_t) * lfdll_set_mark__link_ptr = (is_next) ?                                         \
            LFDLL__LINK_NEXT((lst_ptr), (node_ptr)) :                                                       \
            LFDLL__LINK_PREV((lst_ptr), (node_ptr));                                                        \
                                                                                                            \
        uintptr_t lfdll_set_mark__link1 = atomic_load(lfdll_set_mark__link_ptr);                            \
        while ( ! LFDLL_IS_MARKED(lfdll_set_mark__link1) ) {                                                \
            if(LFDLL__CAS_LINK(lfdll_set_mark__link_ptr,                                                    \
                                lfdll_set_mark__link1,                                                      \
                                LFDLL_MARK(lfdll_set_mark__link1)))                                         \
                break;                                                                                      \
            LFDLL_BACKOFF();                                                                                \
            lfdll_set_mark__link1 = atomic_load(lfdll_set_mark__link_ptr);                                  \
        }                                                                                                   \
    } while( 0 )



/************************************************************************
 *
 * LFDLL__FIND_PRED_FORWARD_BOUNDED()
 *
 *     Perform a bounded forward scan from head_s to find the predecessor of
 *     node_ptr.
 *
 *     This helper is used only as a bounded recovery path when a direct
 *     backward repair attempt has taken too many hops.
 *
 ************************************************************************/
#define LFDLL__FIND_PRED_FORWARD_BOUNDED(lst_ptr, node_ptr, max_steps, pred_out)                            \
    do {                                                                                                    \
        void * lfdll_fpfb__pred = (void *)&((lst_ptr)->head_s);                                             \
        void * lfdll_fpfb__cur  = NULL;                                                                     \
        unsigned int lfdll_fpfb__steps = 0U;                                                                \
        bool lfdll_fpfb__found = false;                                                                     \
                                                                                                            \
        assert(lst_ptr);                                                                                    \
        assert(node_ptr);                                                                                   \
        assert(pred_out);                                                                                   \
                                                                                                            \
        *(pred_out) = (void *)&((lst_ptr)->head_s);                                                         \
                                                                                                            \
        LFDLL__LOAD_PTR((lst_ptr), lfdll_fpfb__pred, true, &lfdll_fpfb__cur);                               \
                                                                                                            \
        while (lfdll_fpfb__steps++ < (unsigned int)(max_steps)) {                                           \
            if(lfdll_fpfb__cur == (node_ptr)) {                                                             \
                *(pred_out) = lfdll_fpfb__pred;                                                             \
                lfdll_fpfb__found = true;                                                                   \
                break;                                                                                      \
            }                                                                                               \
            if(lfdll_fpfb__cur == (void *)&((lst_ptr)->tail_s))                                             \
                break;                                                                                      \
            lfdll_fpfb__pred = lfdll_fpfb__cur;                                                             \
            LFDLL__LOAD_PTR((lst_ptr), lfdll_fpfb__cur, true, &lfdll_fpfb__cur);                            \
        }                                                                                                   \
                                                                                                            \
        if(!lfdll_fpfb__found)                                                                              \
            *(pred_out) = (void *)&((lst_ptr)->head_s);                                                     \
    } while( 0 )

/************************************************************************
 *
 * LFDLL__CORRECT_PREV()
 *
 *     Repair node_ptr->prev by walking the forward chain from prev_ptr.
 *
 *     This macro implements the CorrectPrev helper from the Sundell-Tsigas
 *     algorithm. It advances through the authoritative next links, splices out
 *     marked predecessors when possible, and attempts to publish the repaired
 *     predecessor in node_ptr->prev.
 *
 *     On completion, pred_out receives the best live predecessor located for
 *     node_ptr.
 *
 ************************************************************************/
#define LFDLL__CORRECT_PREV(lst_ptr, prev_ptr, node_ptr, pred_out)                                          \
    do {                                                                                                    \
        void * lfdll_cp__prev = (prev_ptr);                                                                 \
        void * lfdll_cp__node = (node_ptr);                                                                 \
        void * lfdll_cp__lastlink = NULL;                                                                   \
        unsigned int lfdll_cp__spins = 0U;                                                                  \
        unsigned int lfdll_cp__hops  = 0U;                                                                  \
        void * lfdll_cp__ret = NULL;                                                                        \
                                                                                                            \
        assert((lst_ptr) && lfdll_cp__node);                                                                \
        assert(prev_ptr);                                                                                   \
        assert(pred_out);                                                                                   \
                                                                                                            \
        atomic_fetch_add(&((lst_ptr)->correct_prev_calls), 1LL);                                            \
                                                                                                            \
        while(true) {                                                                                       \
            lfdll_cp__spins++;                                                                              \
            if(++lfdll_cp__hops >= LFDLL_CORRECT_PREV_MAX_HOPS) {                                           \
                void * lfdll_cp__fpred = NULL;                                                              \
                bool lfdll_cp__ignored_marked;                                                              \
                uintptr_t lfdll_cp__cur_prev_raw;                                                           \
                LFDLL__FIND_PRED_FORWARD_BOUNDED((lst_ptr), lfdll_cp__node, LFDLL_CORRECT_PREV_MAX_HOPS,    \
                                                                            &lfdll_cp__fpred);              \
                LFDLL__LOAD_LINK((lst_ptr), lfdll_cp__node, false, &lfdll_cp__ignored_marked,               \
                                                                   &lfdll_cp__cur_prev_raw);                \
                (void)LFDLL__CAS_LINK(LFDLL__LINK_PREV((lst_ptr), lfdll_cp__node), lfdll_cp__cur_prev_raw,  \
                                                                        (uintptr_t)lfdll_cp__fpred);        \
                lfdll_cp__ret = lfdll_cp__fpred;                                                            \
                break;                                                                                      \
            }                                                                                               \
            if((lfdll_cp__spins & 0xFFU) == 0U) {                                                           \
                LFDLL_BACKOFF();                                                                            \
            }                                                                                               \
                                                                                                            \
            if(lfdll_cp__prev == (void *)&((lst_ptr)->tail_s) &&                                            \
                                lfdll_cp__node != (void *)&((lst_ptr)->tail_s)) {                           \
                lfdll_cp__ret = (void *)&((lst_ptr)->head_s);                                               \
                break;                                                                                      \
            }                                                                                               \
            if(lfdll_cp__spins >= LFDLL_CORRECT_PREV_MAX_SPINS) {                                           \
                lfdll_cp__prev = (void *)&((lst_ptr)->head_s);                                              \
                lfdll_cp__lastlink = NULL;                                                                  \
                lfdll_cp__spins = 0U;                                                                       \
                continue;                                                                                   \
            }                                                                                               \
            atomic_fetch_add(&((lst_ptr)->correct_prev_iters), 1LL);                                        \
            bool lfdll_cp__ignored_marked;                                                                  \
            uintptr_t lfdll_cp__link1;                                                                      \
            LFDLL__LOAD_LINK((lst_ptr), lfdll_cp__node, false, &lfdll_cp__ignored_marked, &lfdll_cp__link1);\
            if(LFDLL_IS_MARKED(lfdll_cp__link1))                                                            \
                break;                                                                                      \
                                                                                                            \
            bool lfdll_cp__prev2_marked = false;                                                            \
            uintptr_t lfdll_cp__prev2_raw;                                                                  \
            void * lfdll_cp__prev2;                                                                         \
            LFDLL__LOAD_LINK((lst_ptr), lfdll_cp__prev, true, &lfdll_cp__prev2_marked,                      \
                                                              &lfdll_cp__prev2_raw);                        \
            lfdll_cp__prev2 = LFDLL_GET_PTR(lfdll_cp__prev2_raw);                                           \
                                                                                                            \
	        if(lfdll_cp__prev2 == lfdll_cp__prev && lfdll_cp__prev2 != lfdll_cp__node) {                    \
                if(lfdll_cp__prev == (void *)&((lst_ptr)->tail_s)) {                                        \
                    lfdll_cp__ret = lfdll_cp__prev;                                                         \
                    break;                                                                                  \
                }                                                                                           \
                lfdll_cp__prev = (void *)&((lst_ptr)->head_s);                                              \
                lfdll_cp__lastlink = NULL;                                                                  \
                continue;                                                                                   \
            }                                                                                               \
                                                                                                            \
            if(lfdll_cp__prev2_marked) {                                                                    \
                if(lfdll_cp__lastlink != NULL) {                                                            \
                    LFDLL__SET_MARK((lst_ptr), lfdll_cp__prev, false);                                      \
                    {                                                                                       \
                        _Atomic uintptr_t * lfdll_cp__ll_next_ptr = LFDLL__LINK_NEXT((lst_ptr),             \
                                                                    lfdll_cp__lastlink);                    \
                        uintptr_t lfdll_cp__ll_raw = atomic_load_explicit(lfdll_cp__ll_next_ptr,            \
                                                                         memory_order_acquire);             \
                                                                                                            \
                        if(LFDLL_GET_PTR(lfdll_cp__ll_raw) == lfdll_cp__prev) {                             \
                            uintptr_t lfdll_cp__desired = ((uintptr_t)lfdll_cp__prev2) |                    \
                                                          (lfdll_cp__ll_raw & 1ULL);                        \
                            uintptr_t lfdll_cp__exp = lfdll_cp__ll_raw;                                     \
                                                                                                            \
                            (void)atomic_compare_exchange_strong_explicit(                                  \
                                lfdll_cp__ll_next_ptr, &lfdll_cp__exp, lfdll_cp__desired,                   \
                                             memory_order_acq_rel, memory_order_acquire);                   \
                        }                                                                                   \
                    }                                                                                       \
                    lfdll_cp__prev = lfdll_cp__lastlink;                                                    \
                    lfdll_cp__lastlink = NULL;                                                              \
                    continue;                                                                               \
                }                                                                                           \
                lfdll_cp__prev = (void *)&((lst_ptr)->head_s);                                              \
                lfdll_cp__lastlink = NULL;                                                                  \
                continue;                                                                                   \
            }                                                                                               \
                                                                                                            \
            if(lfdll_cp__prev2 != lfdll_cp__node) {                                                         \
                lfdll_cp__lastlink = lfdll_cp__prev;                                                        \
                lfdll_cp__prev = lfdll_cp__prev2;                                                           \
                continue;                                                                                   \
            }                                                                                               \
                                                                                                            \
            if(LFDLL__CAS_LINK(LFDLL__LINK_PREV((lst_ptr), lfdll_cp__node), lfdll_cp__link1,                \
                                                                            (uintptr_t)lfdll_cp__prev)) {   \
                bool lfdll_cp__pp_marked = false;                                                           \
                uintptr_t lfdll_cp__tmp;                                                                    \
                LFDLL__LOAD_LINK((lst_ptr), lfdll_cp__prev, false, &lfdll_cp__pp_marked, &lfdll_cp__tmp);   \
                if(lfdll_cp__pp_marked)                                                                     \
                    continue;                                                                               \
                break;                                                                                      \
            }                                                                                               \
        }                                                                                                   \
                                                                                                            \
        if(!lfdll_cp__ret)                                                                                  \
            lfdll_cp__ret = lfdll_cp__prev;                                                                 \
                                                                                                            \
        *(pred_out) = lfdll_cp__ret;                                                                        \
    } while( 0 )


/************************************************************************
 *
 * LFDLL__PUSH_COMMON()
 *
 *     Complete the insertion helping step after linking node_ptr into the
 *     forward chain.
 *
 *     Once prev->next has been swung to node_ptr, this helper attempts to
 *     update next_ptr->prev to point back to node_ptr, or triggers backward
 *     repair if node_ptr is concurrently being deleted.
 *
 ************************************************************************/
#define LFDLL__PUSH_COMMON(lst_ptr, node_ptr, next_ptr)                                                     \
    do {                                                                                                    \
        while(true) {                                                                                       \
            bool lfdll_pc__ignored_marked;                                                                  \
            uintptr_t lfdll_pc__link1;                                                                      \
            LFDLL__LOAD_LINK((lst_ptr), (next_ptr), false, &lfdll_pc__ignored_marked, &lfdll_pc__link1);    \
                                                                                                            \
            if(LFDLL_IS_MARKED(lfdll_pc__link1))                                                            \
                break;                                                                                      \
            uintptr_t lfdll_pc__node_next;                                                                  \
            LFDLL__LOAD_LINK((lst_ptr), (node_ptr), true, &lfdll_pc__ignored_marked, &lfdll_pc__node_next); \
            if(lfdll_pc__node_next != (uintptr_t)(next_ptr))                                                \
                break;                                                                                      \
                                                                                                            \
            if(LFDLL__CAS_LINK(LFDLL__LINK_PREV((lst_ptr), (next_ptr)), lfdll_pc__link1,                    \
                                                                        (uintptr_t)(node_ptr))) {           \
                bool lfdll_pc__np_marked = false;                                                           \
                uintptr_t lfdll_pc__tmp;                                                                    \
                LFDLL__LOAD_LINK((lst_ptr), (node_ptr), false, &lfdll_pc__np_marked, &lfdll_pc__tmp);       \
                if(lfdll_pc__np_marked) {                                                                   \
                    void * lfdll_pc__ignored = NULL;                                                        \
                    LFDLL__CORRECT_PREV((lst_ptr), (node_ptr), (next_ptr), &lfdll_pc__ignored);             \
                }                                                                                           \
                break;                                                                                      \
            }                                                                                               \
        }                                                                                                   \
    } while( 0 )

/************************************************************************
 *
 * Public API macros.
 *
 *     The list exposes a small typeless interface: initialize/create/destroy
 *     the list object, prepend/append/remove user nodes, query the front/back
 *     entries, and dump internal statistics for debugging.
 *
 ************************************************************************/


/************************************************************************
 *
 * LFDLL_INIT()
 *
 *     Initialize an instance of lfdll_list_t.
 *
 *     The caller supplies the byte offsets of the embedded next and prev link
 *     fields in the actual node type. The macro initializes both sentinels,
 *     the associated free-list bookkeeping, and the basic list statistics.
 *
 *                                             CMM -- 3/4/26
 *
 ************************************************************************/
#define LFDLL_INIT(lst_ptr, next_off, prev_off)                                                             \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        void * lfdll_init__lst_vp = (void *)(lst_ptr);                                                      \
                                                                                                            \
        if(!lfdll_init__lst_vp)                                                                             \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
                                                                                                            \
        if((next_off) == (prev_off))                                                                        \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "next_off == prev_off");                              \
                                                                                                            \
        assert(((next_off) % _Alignof(_Atomic(uintptr_t))) == 0U);                                          \
        assert(((prev_off) % _Alignof(_Atomic(uintptr_t))) == 0U);                                          \
                                                                                                            \
        (lst_ptr)->tag      = LFDLL_VALID;                                                                  \
        (lst_ptr)->off_next = (next_off);                                                                   \
        (lst_ptr)->off_prev = (prev_off);                                                                   \
                                                                                                            \
        atomic_store(&((lst_ptr)->head_s.prev), (uintptr_t)&((lst_ptr)->head_s));                           \
        atomic_store(&((lst_ptr)->tail_s.next), (uintptr_t)&((lst_ptr)->tail_s));                           \
        atomic_store(&((lst_ptr)->head_s.next), (uintptr_t)&((lst_ptr)->tail_s));                           \
        atomic_store(&((lst_ptr)->tail_s.prev), (uintptr_t)&((lst_ptr)->head_s));                           \
                                                                                                            \
        LFFL_INIT(&((lst_ptr)->fl));                                                                        \
                                                                                                            \
        atomic_store(&((lst_ptr)->ins_ok), 0LL);                                                            \
        atomic_store(&((lst_ptr)->ins_fail), 0LL);                                                          \
        atomic_store(&((lst_ptr)->rem_ok), 0LL);                                                            \
        atomic_store(&((lst_ptr)->rem_fail), 0LL);                                                          \
        atomic_store(&((lst_ptr)->correct_prev_calls), 0LL);                                                \
        atomic_store(&((lst_ptr)->correct_prev_iters), 0LL);                                                \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
        assert((uintptr_t)&((lst_ptr)->head_s) == atomic_load(&((lst_ptr)->tail_s.prev)));                  \
        assert((uintptr_t)&((lst_ptr)->tail_s) == atomic_load(&((lst_ptr)->head_s.next)));                  \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_CREATE()
 *
 *     Allocate and initialize a new lfdll_list_t on the heap.
 *
 *     LFDLL_INIT() initializes a caller-provided list object. LFDLL_CREATE()
 *     allocates the list object and then invokes LFDLL_INIT() on it.
 *
 *                                             CMM -- 3/4/26
 *
 ************************************************************************/
#define LFDLL_CREATE(next_off, prev_off, out_lst_ptr)                                                       \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        lfdll_list_t * lfdll_create__lst_ptr = NULL;                                                        \
        void * lfdll_create__out_vp = (void *)(out_lst_ptr);                                                \
                                                                                                            \
        if(!lfdll_create__out_vp)                                                                           \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "out_lst_ptr == NULL");                               \
                                                                                                            \
        lfdll_create__lst_ptr = (lfdll_list_t *)LFDLL_MALLOC(sizeof(*lfdll_create__lst_ptr));               \
        if(!lfdll_create__lst_ptr)                                                                          \
            LFDLL_ERROR(H5E_RESOURCE, H5E_CANTALLOC, FAIL, "malloc failed");                                \
                                                                                                            \
        LFDLL_INIT(lfdll_create__lst_ptr, (next_off), (prev_off));                                          \
        assert(LFDLL_VALID == lfdll_create__lst_ptr->tag);                                                  \
        *(out_lst_ptr) = lfdll_create__lst_ptr;                                                             \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_DESTROY()
 *
 *     Free a heap-allocated lfdll_list_t.
 *
 *     This is a shallow destroy of the list object itself. Ownership and
 *     reclamation of user nodes remain the caller's responsibility.
 *
 *                                             CMM -- 3/4/26
 *
 ************************************************************************/
#define LFDLL_DESTROY(lst_ptr)                                                                              \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        void * lfdll_destroy__lst_vp = (void *)(lst_ptr);                                                   \
                                                                                                            \
        if(!lfdll_destroy__lst_vp)                                                                          \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
        assert((uintptr_t)&((lst_ptr)->head_s) == atomic_load(&((lst_ptr)->tail_s.prev)) ||                 \
               (uintptr_t)&((lst_ptr)->head_s) != 0U);                                                      \
                                                                                                            \
        LFDLL_FINI(lst_ptr);                                                                                \
        LFDLL_FREE(lst_ptr);                                                                                \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_PREPEND()
 *
 *     Insert node_ptr at the front of the list, immediately after head_s.
 *
 *     This macro follows the PushLeft-style insertion from the reference
 *     algorithm. It publishes node_ptr through the authoritative forward link,
 *     then completes the backward-link repair/helping step with
 *     LFDLL__PUSH_COMMON().
 *
 *     The operation result is returned through out_rc.
 *
 *                                             CMM -- 3/4/26
 *
 ************************************************************************/
#define LFDLL_PREPEND(lst_ptr, node_ptr, out_rc)                                                            \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        int lfdll_prepend__rc = FAIL;                                                                       \
        int * lfdll_prepend__out_rc_ptr = (out_rc);                                                         \
        void * lfdll_prepend__out_vp = (void *)lfdll_prepend__out_rc_ptr;                                   \
        void * lfdll_prepend__lst_vp = (void *)(lst_ptr);                                                   \
        void * lfdll_prepend__node_vp = (void *)(node_ptr);                                                 \
                                                                                                            \
        if(!lfdll_prepend__lst_vp)                                                                          \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
        if(!lfdll_prepend__node_vp)                                                                         \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "node_ptr == NULL");                                  \
        if(!lfdll_prepend__out_vp)                                                                          \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "out_rc == NULL");                                    \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
        assert((node_ptr) != (void *)&((lst_ptr)->head_s));                                                 \
        assert((node_ptr) != (void *)&((lst_ptr)->tail_s));                                                 \
        assert((((uintptr_t)(node_ptr)) & LFDLL_TAG_MASK) == 0U);                                           \
                                                                                                            \
        lffl_token_t * lfdll_prepend__tok_ptr = NULL;                                                       \
        LFFL_ENTER(&((lst_ptr)->fl), &lfdll_prepend__tok_ptr);                                              \
                                                                                                            \
        void * lfdll_prepend__prev = (void *)&((lst_ptr)->head_s);                                          \
                                                                                                            \
        while(true) {                                                                                       \
            void * lfdll_prepend__next = NULL;                                                              \
            LFDLL__LOAD_PTR((lst_ptr), lfdll_prepend__prev, true, &lfdll_prepend__next);                    \
                                                                                                            \
            atomic_store(LFDLL__LINK_PREV((lst_ptr), (node_ptr)), (uintptr_t)lfdll_prepend__prev);          \
            atomic_store(LFDLL__LINK_NEXT((lst_ptr), (node_ptr)), (uintptr_t)lfdll_prepend__next);          \
                                                                                                            \
            if(LFDLL__CAS_LINK(LFDLL__LINK_NEXT((lst_ptr), lfdll_prepend__prev),                            \
                                (uintptr_t)lfdll_prepend__next,                                             \
                                (uintptr_t)(node_ptr))) {                                                   \
                LFDLL__PUSH_COMMON((lst_ptr), (node_ptr), lfdll_prepend__next);                             \
                LFFL_EXIT(&((lst_ptr)->fl), lfdll_prepend__tok_ptr);                                        \
                atomic_fetch_add(&((lst_ptr)->ins_ok), 1LL);                                                \
                lfdll_prepend__rc = SUCCEED;                                                                \
                break;                                                                                      \
            }                                                                                               \
            LFDLL_BACKOFF();                                                                                \
        }                                                                                                   \
                                                                                                            \
        *(lfdll_prepend__out_rc_ptr) = lfdll_prepend__rc;                                                   \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_APPEND()
 *
 *     Insert node_ptr at the back of the list, immediately before tail_s.
 *
 *     This macro follows the PushRight-style insertion from the reference
 *     algorithm. The predecessor of tail_s is located through the forward
 *     chain, node_ptr is linked through prev->next, and the backward-link help
 *     step is completed with LFDLL__PUSH_COMMON().
 *
 *     The operation result is returned through out_rc.
 *
 *                                             CMM -- 3/3/26
 *
 ************************************************************************/
#define LFDLL_APPEND(lst_ptr, node_ptr, out_rc)                                                             \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        int lfdll_append__rc = FAIL;                                                                        \
        int * lfdll_append__out_rc_ptr = (out_rc);                                                          \
        void * lfdll_append__out_vp = (void *)lfdll_append__out_rc_ptr;                                     \
        void * lfdll_append__lst_vp = (void *)(lst_ptr);                                                    \
        void * lfdll_append__node_vp = (void *)(node_ptr);                                                  \
                                                                                                            \
        if(!lfdll_append__lst_vp)                                                                           \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
        if(!lfdll_append__node_vp)                                                                          \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "node_ptr == NULL");                                  \
        if(!lfdll_append__out_vp)                                                                           \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "out_rc == NULL");                                    \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
        assert((node_ptr) != (void *)&((lst_ptr)->head_s));                                                 \
        assert((node_ptr) != (void *)&((lst_ptr)->tail_s));                                                 \
        assert((((uintptr_t)(node_ptr)) & LFDLL_TAG_MASK) == 0U);                                           \
                                                                                                            \
        lffl_token_t * lfdll_append__tok_ptr = NULL;                                                        \
        LFFL_ENTER(&((lst_ptr)->fl), &lfdll_append__tok_ptr);                                               \
                                                                                                            \
        void * lfdll_append__next = (void *)&((lst_ptr)->tail_s);                                           \
        void * lfdll_append__prev = NULL;                                                                   \
                                                                                                            \
        while(true) {                                                                                       \
            LFDLL__FIND_PRED_FORWARD_BOUNDED((lst_ptr), lfdll_append__next,                                 \
                                            LFDLL_CORRECT_PREV_MAX_HOPS, &lfdll_append__prev);              \
                                                                                                            \
            bool lfdll_append__pn_marked = false;                                                           \
            uintptr_t lfdll_append__pn_raw;                                                                 \
            LFDLL__LOAD_LINK((lst_ptr), lfdll_append__prev, true, &lfdll_append__pn_marked,                 \
                                                                  &lfdll_append__pn_raw);                   \
            if(lfdll_append__pn_marked) {                                                                   \
                void * lfdll_append__tmp = NULL;                                                            \
                LFDLL__CORRECT_PREV((lst_ptr), lfdll_append__prev, lfdll_append__next, &lfdll_append__tmp); \
                continue;                                                                                   \
            }                                                                                               \
                                                                                                            \
            atomic_store(LFDLL__LINK_PREV((lst_ptr), (node_ptr)), (uintptr_t)lfdll_append__prev);           \
            atomic_store(LFDLL__LINK_NEXT((lst_ptr), (node_ptr)), (uintptr_t)lfdll_append__next);           \
                                                                                                            \
            if(LFDLL__CAS_LINK(LFDLL__LINK_NEXT((lst_ptr), lfdll_append__prev),                             \
                                (uintptr_t)lfdll_append__next,                                              \
                                (uintptr_t)(node_ptr))) {                                                   \
                LFDLL__PUSH_COMMON((lst_ptr), (node_ptr), lfdll_append__next);                              \
                LFFL_EXIT(&((lst_ptr)->fl), lfdll_append__tok_ptr);                                         \
                atomic_fetch_add(&((lst_ptr)->ins_ok), 1LL);                                                \
                lfdll_append__rc = SUCCEED;                                                                 \
                break;                                                                                      \
            }                                                                                               \
            LFDLL_BACKOFF();                                                                                \
        }                                                                                                   \
                                                                                                            \
        *(lfdll_append__out_rc_ptr) = lfdll_append__rc;                                                     \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_REMOVE()
 *
 *     Remove node_ptr from the list by pointer.
 *
 *     This macro follows the Delete operation from the reference algorithm.
 *     It first marks node_ptr->next to logically delete the node, then marks
 *     node_ptr->prev, helps repair neighboring backward links, and optionally
 *     retires the node through the free-list helper when deleter_fn is not
 *     NULL.
 *
 *     The operation result is returned through out_rc.
 *
 *                                             CMM -- 3/3/26
 *
 ************************************************************************/
#define LFDLL_REMOVE(lst_ptr, node_ptr, deleter_fn, out_rc)                                                 \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        int lfdll_remove__rc = FAIL;                                                                        \
        int * lfdll_remove__out_rc_ptr = (out_rc);                                                          \
        void * lfdll_remove__out_vp = (void *)lfdll_remove__out_rc_ptr;                                     \
        void * lfdll_remove__lst_vp = (void *)(lst_ptr);                                                    \
        void * lfdll_remove__node_vp = (void *)(node_ptr);                                                  \
        bool lfdll_remove__early_fail = false;                                                              \
                                                                                                            \
        if(!lfdll_remove__lst_vp)                                                                           \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
        if(!lfdll_remove__node_vp)                                                                          \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "node_ptr == NULL");                                  \
        if(!lfdll_remove__out_vp)                                                                           \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "out_rc == NULL");                                    \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
        assert((((uintptr_t)(node_ptr)) & LFDLL_TAG_MASK) == 0U);                                           \
                                                                                                            \
        if((node_ptr) == (void *)&((lst_ptr)->head_s) || (node_ptr) == (void *)&((lst_ptr)->tail_s)) {      \
            *(lfdll_remove__out_rc_ptr) = FAIL;                                                             \
            break;                                                                                          \
        }                                                                                                   \
                                                                                                            \
        lffl_token_t * lfdll_remove__tok_ptr = NULL;                                                        \
        LFFL_ENTER(&((lst_ptr)->fl), &lfdll_remove__tok_ptr);                                               \
                                                                                                            \
        while(true) {                                                                                       \
            bool lfdll_remove__ignored_marked;                                                              \
            uintptr_t lfdll_remove__next_raw;                                                               \
            LFDLL__LOAD_LINK((lst_ptr), (node_ptr), true, &lfdll_remove__ignored_marked,                    \
                                                          &lfdll_remove__next_raw);                         \
                                                                                                            \
            if(LFDLL_IS_MARKED(lfdll_remove__next_raw)) {                                                   \
                LFFL_EXIT(&((lst_ptr)->fl), lfdll_remove__tok_ptr);                                         \
                *(lfdll_remove__out_rc_ptr) = FAIL;                                                         \
                lfdll_remove__early_fail = true;                                                            \
                break;                                                                                      \
            }                                                                                               \
                                                                                                            \
            if(LFDLL__CAS_LINK(LFDLL__LINK_NEXT((lst_ptr), (node_ptr)),                                     \
                                lfdll_remove__next_raw,                                                     \
                                LFDLL_MARK(lfdll_remove__next_raw)))                                        \
                break;                                                                                      \
            LFDLL_BACKOFF();                                                                                \
        }                                                                                                   \
                                                                                                            \
        if(lfdll_remove__early_fail)                                                                        \
            break;                                                                                          \
                                                                                                            \
        while(true) {                                                                                       \
            bool lfdll_remove__ignored_marked;                                                              \
            uintptr_t lfdll_remove__prev_raw;                                                               \
            LFDLL__LOAD_LINK((lst_ptr), (node_ptr), false, &lfdll_remove__ignored_marked,                   \
                                                           &lfdll_remove__prev_raw);                        \
                                                                                                            \
            if(LFDLL_IS_MARKED(lfdll_remove__prev_raw))                                                     \
                break;                                                                                      \
            if(LFDLL__CAS_LINK(LFDLL__LINK_PREV((lst_ptr), (node_ptr)),                                     \
                                lfdll_remove__prev_raw,                                                     \
                                LFDLL_MARK(lfdll_remove__prev_raw)))                                        \
                break;                                                                                      \
            LFDLL_BACKOFF();                                                                                \
        }                                                                                                   \
                                                                                                            \
        void * lfdll_remove__prev = NULL;                                                                   \
        void * lfdll_remove__next = NULL;                                                                   \
        LFDLL__LOAD_PTR((lst_ptr), (node_ptr), false, &lfdll_remove__prev);                                 \
        LFDLL__LOAD_PTR((lst_ptr), (node_ptr), true, &lfdll_remove__next);                                  \
                                                                                                            \
        unsigned int lfdll_remove__skip = 0U;                                                               \
        while (lfdll_remove__next != (void *)&((lst_ptr)->tail_s)) {                                        \
            if(++lfdll_remove__skip >= 1024U) {                                                             \
                lfdll_remove__next = (void *)&((lst_ptr)->tail_s);                                          \
                break;                                                                                      \
            }                                                                                               \
            bool lfdll_remove__nprev_marked = false;                                                        \
            uintptr_t lfdll_remove__tmp_raw;                                                                \
            LFDLL__LOAD_LINK((lst_ptr), lfdll_remove__next, false, &lfdll_remove__nprev_marked,             \
                                                                   &lfdll_remove__tmp_raw);                 \
            if(!lfdll_remove__nprev_marked)                                                                 \
                break;                                                                                      \
            void * lfdll_remove__nnext = NULL;                                                              \
            LFDLL__LOAD_PTR((lst_ptr), lfdll_remove__next, true, &lfdll_remove__nnext);                     \
            if(lfdll_remove__nnext == lfdll_remove__next) {                                                 \
                lfdll_remove__next = (void *)&((lst_ptr)->tail_s);                                          \
                break;                                                                                      \
            }                                                                                               \
            lfdll_remove__next = lfdll_remove__nnext;                                                       \
        }                                                                                                   \
                                                                                                            \
        if(lfdll_remove__next != (void *)&((lst_ptr)->tail_s)) {                                            \
            bool lfdll_remove__ignored_marked;                                                              \
            uintptr_t lfdll_remove__cur_prev_raw;                                                           \
            LFDLL__LOAD_LINK((lst_ptr), lfdll_remove__next, false, &lfdll_remove__ignored_marked,           \
                                                                   &lfdll_remove__cur_prev_raw);            \
            if(LFDLL_PTR(lfdll_remove__cur_prev_raw) == (node_ptr) &&                                       \
                                                        !LFDLL_IS_MARKED(lfdll_remove__cur_prev_raw)) {     \
                uintptr_t lfdll_remove__cur_next_raw;                                                       \
                LFDLL__LOAD_LINK((lst_ptr), lfdll_remove__prev, true, &lfdll_remove__ignored_marked,        \
                                                                      &lfdll_remove__cur_next_raw);         \
                if(LFDLL_PTR(lfdll_remove__cur_next_raw) == lfdll_remove__next &&                           \
                                                            !LFDLL_IS_MARKED(lfdll_remove__cur_next_raw))   \
                    (void)LFDLL__CAS_LINK(LFDLL__LINK_PREV((lst_ptr), lfdll_remove__next),                  \
                                         lfdll_remove__cur_prev_raw,                                        \
                                         (uintptr_t)lfdll_remove__prev);                                    \
            }                                                                                               \
        }                                                                                                   \
                                                                                                            \
        {                                                                                                   \
            void * lfdll_remove__ignored = NULL;                                                            \
            LFDLL__CORRECT_PREV((lst_ptr), lfdll_remove__prev, lfdll_remove__next, &lfdll_remove__ignored); \
        }                                                                                                   \
                                                                                                            \
        {                                                                                                   \
            unsigned int lfdll_remove__rounds = 0U;                                                         \
            while (lfdll_remove__rounds++ < LFDLL_REMOVE_WAIT_MAX_ROUNDS) {                                 \
                void * lfdll_remove__pred = (void *)&((lst_ptr)->head_s);                                   \
                void * lfdll_remove__cur = NULL;                                                            \
                unsigned int lfdll_remove__steps = 0U;                                                      \
                bool lfdll_remove__found = false;                                                           \
                                                                                                            \
                LFDLL__LOAD_PTR((lst_ptr), lfdll_remove__pred, true, &lfdll_remove__cur);                   \
                                                                                                            \
                while (lfdll_remove__cur != (void *)&((lst_ptr)->tail_s) &&                                 \
                       lfdll_remove__steps++ < LFDLL_REMOVE_WAIT_SCAN_LIMIT) {                              \
                    lfdll_remove__cur = LFDLL_PTR((uintptr_t)lfdll_remove__cur);                            \
                    if(lfdll_remove__cur == (node_ptr)) {                                                   \
                        lfdll_remove__found = true;                                                         \
                        break;                                                                              \
                    }                                                                                       \
                    lfdll_remove__pred = lfdll_remove__cur;                                                 \
                    LFDLL__LOAD_PTR((lst_ptr), lfdll_remove__cur, true, &lfdll_remove__cur);                \
                }                                                                                           \
                                                                                                            \
                if(!lfdll_remove__found)                                                                    \
                    break;                                                                                  \
                LFDLL_BACKOFF();                                                                            \
                                                                                                            \
                void * lfdll_remove__succ = NULL;                                                           \
                LFDLL__LOAD_PTR((lst_ptr), (node_ptr), true, &lfdll_remove__succ);                          \
                lfdll_remove__succ = LFDLL_PTR((uintptr_t)lfdll_remove__succ);                              \
                {                                                                                           \
                    void * lfdll_remove__tmp = NULL;                                                        \
                    LFDLL__CORRECT_PREV((lst_ptr), lfdll_remove__pred, lfdll_remove__succ,                  \
                                                                       &lfdll_remove__tmp);                 \
                }                                                                                           \
            }                                                                                               \
        }                                                                                                   \
                                                                                                            \
        if(deleter_fn)                                                                                      \
            LFFL_RETIRE(&((lst_ptr)->fl), (node_ptr), (deleter_fn));                                        \
                                                                                                            \
        LFFL_EXIT(&((lst_ptr)->fl), lfdll_remove__tok_ptr);                                                 \
        atomic_fetch_add(&((lst_ptr)->rem_ok), 1LL);                                                        \
        lfdll_remove__rc = SUCCEED;                                                                         \
        *(lfdll_remove__out_rc_ptr) = lfdll_remove__rc;                                                     \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_FRONT()
 *
 *     Return the first live user node in the list, or NULL if the list is
 *     empty.
 *
 *     The traversal skips over nodes whose forward links are already marked
 *     for logical deletion.
 *
 *                                             CMM -- 3/3/26
 *
 ************************************************************************/
#define LFDLL_FRONT(lst_ptr, out_node_ptr)                                                                  \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        void ** lfdll_front__out_node_ptr = (void **)(out_node_ptr);                                        \
        void * lfdll_front__out_vp = (void *)lfdll_front__out_node_ptr;                                     \
        void * lfdll_front__lst_vp = (void *)(lst_ptr);                                                     \
        void * lfdll_front__res = NULL;                                                                     \
                                                                                                            \
        if(!lfdll_front__lst_vp)                                                                            \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
        if(!lfdll_front__out_vp)                                                                            \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "out_node_ptr == NULL");                              \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
                                                                                                            \
        void * lfdll_front__node = NULL;                                                                    \
        LFDLL__LOAD_PTR((lst_ptr), &((lst_ptr)->head_s), true, &lfdll_front__node);                         \
                                                                                                            \
        while (lfdll_front__node != (void *)&((lst_ptr)->tail_s)) {                                         \
            bool lfdll_front__marked = false;                                                               \
            uintptr_t lfdll_front__tmp;                                                                     \
            LFDLL__LOAD_LINK((lst_ptr), lfdll_front__node, true, &lfdll_front__marked, &lfdll_front__tmp);  \
            if(!lfdll_front__marked) {                                                                      \
                lfdll_front__res = lfdll_front__node;                                                       \
                break;                                                                                      \
            }                                                                                               \
            LFDLL__LOAD_PTR((lst_ptr), lfdll_front__node, true, &lfdll_front__node);                        \
        }                                                                                                   \
                                                                                                            \
        *(lfdll_front__out_node_ptr) = lfdll_front__res;                                                    \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_BACK()
 *
 *     Return the last live user node in the list, or NULL if the list is
 *     empty.
 *
 *     Because backward links are advisory, the macro first repairs the tail
 *     predecessor path with LFDLL__CORRECT_PREV() before reporting the result.
 *
 *                                             CMM -- 3/4/26
 *
 ************************************************************************/
#define LFDLL_BACK(lst_ptr, out_node_ptr)                                                                   \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        void ** lfdll_back__out_node_ptr = (void **)(out_node_ptr);                                         \
        void * lfdll_back__out_vp = (void *)lfdll_back__out_node_ptr;                                       \
        void * lfdll_back__lst_vp = (void *)(lst_ptr);                                                      \
        void * lfdll_back__res = NULL;                                                                      \
                                                                                                            \
        if(!lfdll_back__lst_vp)                                                                             \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
        if(!lfdll_back__out_vp)                                                                             \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "out_node_ptr == NULL");                              \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
                                                                                                            \
        void * lfdll_back__node = NULL;                                                                     \
        LFDLL__LOAD_PTR((lst_ptr), &((lst_ptr)->tail_s), false, &lfdll_back__node);                         \
        if(lfdll_back__node == (void *)&((lst_ptr)->head_s)) {                                              \
            *(lfdll_back__out_node_ptr) = NULL;                                                             \
            break;                                                                                          \
        }                                                                                                   \
                                                                                                            \
        void * lfdll_back__tmp = NULL;                                                                      \
        LFDLL__CORRECT_PREV((lst_ptr), lfdll_back__node, &((lst_ptr)->tail_s), &lfdll_back__tmp);           \
                                                                                                            \
        void * lfdll_back__pred = NULL;                                                                     \
        LFDLL__LOAD_PTR((lst_ptr), &((lst_ptr)->tail_s), false, &lfdll_back__pred);                         \
                                                                                                            \
        if(lfdll_back__pred != (void *)&((lst_ptr)->head_s)) {                                              \
            bool lfdll_back__marked = false;                                                                \
            uintptr_t lfdll_back__tmp_raw;                                                                  \
            LFDLL__LOAD_LINK((lst_ptr), lfdll_back__pred, true, &lfdll_back__marked, &lfdll_back__tmp_raw); \
            if(!lfdll_back__marked)                                                                         \
                lfdll_back__res = lfdll_back__pred;                                                         \
        }                                                                                                   \
                                                                                                            \
        *(lfdll_back__out_node_ptr) = lfdll_back__res;                                                      \
    } while( 0 )                                                                                            \



/************************************************************************
 *
 * LFDLL_DUMP_STATS()
 *
 *     Dump basic list statistics to file_ptr.
 *
 *     This is a debugging/reporting helper only. It does not attempt to
 *     serialize concurrent updates to the counters beyond the atomic loads.
 *
 *                                             CMM -- 3/4/26
 *
 ************************************************************************/
#define LFDLL_DUMP_STATS(lst_ptr, file_ptr)                                                                 \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        void * lfdll_dump__lst_vp = (void *)(lst_ptr);                                                      \
        FILE * lfdll_dump__file_ptr = (file_ptr);                                                           \
        void * lfdll_dump__file_vp = (void *)lfdll_dump__file_ptr;                                          \
                                                                                                            \
        if(!lfdll_dump__lst_vp)                                                                             \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
        if(!lfdll_dump__file_vp)                                                                            \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "file_ptr == NULL");                                  \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
                                                                                                            \
        fprintf((file_ptr),                                                                                 \
                "LFDLL stats: ins_ok=%lld ins_fail=%lld rem_ok=%lld rem_fail=%lld cp_calls=%lld             \
                                                                                  cp_iters=%lld\n",         \
                (long long)atomic_load(&((lst_ptr)->ins_ok)),                                               \
                (long long)atomic_load(&((lst_ptr)->ins_fail)),                                             \
                (long long)atomic_load(&((lst_ptr)->rem_ok)),                                               \
                (long long)atomic_load(&((lst_ptr)->rem_fail)),                                             \
                (long long)atomic_load(&((lst_ptr)->correct_prev_calls)),                                   \
                (long long)atomic_load(&((lst_ptr)->correct_prev_iters)));                                  \
    } while( 0 )


/************************************************************************
 *
 * LFDLL_FINI()
 *
 *     Finalize a caller-provided list object without freeing the list
 *     structure itself.
 * 
 *     This macro is intended for lists whose storage is owned by the
 *     caller, such as stack-allocated instances or list objects embedded
 *     in a larger structure. It releases internal LFDLL resources, tears
 *     down the associated free-list state, and marks the list object as
 *     no longer valid for use. The outer lfdll_t object itself is not
 *     freed.
 *
 *     After LFDLL_FINI() returns, the caller must not perform further
 *     LFDLL operations on the list unless the list is initialized again.
 *
 ************************************************************************/
#define LFDLL_FINI(lst_ptr)                                                                                 \
    do {                                                                                                    \
        int ret_value = 0;                                                                                  \
        (void)ret_value;                                                                                    \
        void * lfdll_fini__lst_vp = (void *)(lst_ptr);                                                      \
                                                                                                            \
        if(!lfdll_fini__lst_vp)                                                                             \
            LFDLL_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "lst_ptr == NULL");                                   \
                                                                                                            \
        assert(LFDLL_VALID == (lst_ptr)->tag);                                                              \
                                                                                                            \
        LFFL_DESTROY(&((lst_ptr)->fl));                                                                     \
                                                                                                            \
        (lst_ptr)->tag = 0U;                                                                                \
    } while( 0 ) /* LFDLL_FINI() */

#endif /* LFDLL_H */

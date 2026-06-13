#define _POSIX_C_SOURCE 199309L /* Defines CLOCK_REALTIME */

#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>
#include "lfdll.h"
#include "dispatch_test.h"

/* Macroized usage message for readibility */
#define USAGE_MSG printf("usage: test_dispatch_no_cache [-t <num_threads>] [-s <seed>] [-n <num_requests>] [-v]\n  - num_threads must be greater than 0\n  - num_requests must be greater than 0\n")

/* The dispatch LFDLL will not actually contain request data. */
/* Instead, its purpose is to allow processes to sleep until the element it pushed reaches the front of the queue. */
/* When an element is at the front, the process that owns it wakes and fulfills the request. */
typedef struct queue_element_t {
    struct queue_element_t *next;
    struct queue_element_t *prev;
} queue_element_t;

/* An I/O request will be simulated with a type (0 for read, 1 for write) and a duration. */
typedef struct {
    int type;
    float duration;
} request;

/* Each thread must be initialized with enough data to sustain itself until the end of operation. */
/* num_requests, num_threads, and verbose are all user-defined so they are passed instead of defined globally. */
typedef struct {
    int id;
    pthread_mutex_t *lock;
    pthread_cond_t *cond;
    pthread_mutex_t *write_lock;
    pthread_cond_t *write_cond;
    _Atomic int *write;
    lfdll_list_t *queue;
    request *requests;
    int num_requests;
    int num_threads;
    int verbose;
} thread_info;

long int get_ms_since_epoch() {
    struct timespec ts;

    if (clock_gettime(CLOCK_REALTIME, &ts) == -1)
        return -1;

    return (long int)ts.tv_sec * 1000 + (long int)ts.tv_nsec / 1000000;
}

void *thread_worker(void *vinfo) {
    /* Passing the thread info as a void* and then casting it is necessary for error handling on thread join. */
    thread_info *info = (thread_info *)vinfo;

    /* Sanity check that info was populated. id and verbose are asserted by refrence since they can have a value of 0. */
    assert(info);
    assert(&info->id);
    assert(info->lock);
    assert(info->cond);
    assert(info->write_lock);
    assert(info->write_lock);
    assert(info->write);
    assert(info->queue);
    assert(info->requests);
    assert(info->num_requests);
    assert(info->num_threads);
    assert(&info->verbose);

    /* The thread id (0 <= id < num_threads) will act as an offset for iterating over the requests array. */
    /* A thread will work on all requests with (index % num_threads = id). */
    int i = info->id;

    /* As above, i % num_threads = id. Every request with index (num_threads * n) + id is processed as n increases. */
    while (i < info->num_requests) {
        if (info->verbose)
            printf("Processing request element %d:%d\n", info->id, i);

        /* Append placeholder request to LFDLL. Every thread will have at most one request queued at a time, */
        /*     and will wait for that request to appear at the front of the queue before continuing. */
        queue_element_t *req = calloc(1, sizeof(queue_element_t));
        int *append_result = calloc(1, sizeof(int));
        LFDLL_APPEND(info->queue, req, append_result);
        if (*append_result < 0) {
            if (info->verbose)
                printf("LFDLL_APPEND_ERR: Unable to append.");
            return NULL;
        }
        free(append_result);

        /* Waiting for a condition variable requires the lock be obtained beforehand, so blocking here. */
        pthread_mutex_lock(info->lock);

        /* If the request placed on the LFDLL is at the front, continue, else wait for another signal broadcast and check again. */
        while (true) {
            queue_element_t *front_result = calloc(1, sizeof(queue_element_t));
            LFDLL_FRONT(info->queue, front_result);
            if (front_result->next == req)
                break;
            if (info->verbose)
                printf("Thread %d:%d waiting for front of queue.\n", info->id, i);

            pthread_cond_wait(info->cond, info->lock);
            free(front_result);
        }

        if (info->verbose)
            printf("Removing element %d:%d from queue\n", info->id, i);

        /* There is no pop operation in the LFDLL, so a remove is performed instead */
        /*     with the knowledge that req is at the front of the LFDLL. */
        int *remove_result = calloc(1, sizeof(int));
        LFDLL_REMOVE(info->queue, req, NULL, remove_result);
        if (*remove_result < 0) 
            return NULL;
        free(remove_result);

        /* The front of the LFDLL has changed, signal waiting threads to see if they are next. */
        pthread_mutex_unlock(info->lock);
        pthread_cond_broadcast(info->cond);

        if (info->verbose)
            printf("Thread %d:%d at front of queue.\n", info->id, i);

        /* Simulating I/O with SWMR. */
        /* This request is a write, so no other writes may be active when it is released. */
        if (info->requests[i].type == 1) {
            pthread_mutex_lock(info->write_lock);
            while(true) {
                if (atomic_load(info->write) == 0) { /* Wait until there are no writes happening. */
                    break;
                }

                pthread_cond_wait(info->write_cond, info->write_lock);
            }
            pthread_mutex_unlock(info->write_lock);

            /* Perform write operation. */
            atomic_fetch_add(info->write, 1);
            #ifdef _WIN32
                Sleep((int) info->requests[i].duration * 1000.0);
            #else
                struct timespec sleep = {0, info->requests[i].duration * 1000000000};
                nanosleep(&sleep, NULL);
            #endif
            /* The operation has terminated, decrease the count and signal waiting writes. */
            assert(atomic_fetch_sub(info->write, 1) > 0); /* If we ever try to subtract from 0, error out. */
            pthread_cond_broadcast(info->write_cond);
        } else {
            #ifdef _WIN32
                Sleep((int) info->requests[i].duration * 1000.0);
            #else
                struct timespec sleep = {0, info->requests[i].duration * 1000000000};
                nanosleep(&sleep, NULL);
            #endif
        }

        /* Following the above sleep, the I/O request has successfully been simulated. */
        i += info->num_threads;
    }

    if (info->verbose)
        printf("Thread %d exited.\n", info->id);

    return (void *)1;
}

int main(int argc, char *argv[]) {

    /* Initialize user-defined values to a default */
    long int seed = get_ms_since_epoch();
    int num_threads = 32;
    int num_requests = 1000;
    int verbose = 0;

    /* Parse arguments, looking for flags and their values if necessary. */
    /*     If parsing fails, output usage message and bail. */
    for (int a=1; a<argc; a++) {
        if (argv[a][0] == '-') {
            if (argv[a][1] == 't') {
                if (argv[a][2] != 0) {
                    USAGE_MSG;
                    return 1;
                }
                if (argv[a+1] == NULL) {
                    USAGE_MSG;
                    return 1;
                }

                num_threads = atoi(argv[a+1]);
                if (num_threads < 1 ) {
                    USAGE_MSG;
                    return 1;
                }

                a++;
            } else if (argv[a][1] == 's') {
                if (argv[a][2] != 0) {
                    USAGE_MSG;
                    return 1;
                }
                if (argv[a+1] == NULL) {
                    USAGE_MSG;
                    return 1;
                }

                char *endptr;
                seed = strtol(argv[a+1], &endptr, 10);

                a++;

            } else if (argv[a][1] == 'v') {
                verbose = 1;
            } else if (argv[a][1] == 'n') {
                if (argv[a][2] != 0) {
                    USAGE_MSG;
                    return 1;
                }
                if (argv[a+1] == NULL) {
                    USAGE_MSG;
                    return 1;
                }

                num_requests = atoi(argv[a+1]);

                if (num_requests < 1) {
                    USAGE_MSG;
                    return 1;
                }

                a++;
            }
        } else {
            USAGE_MSG;
            return 1;
        }
    }

    /* Seed the random number generator with seed. This allows reproducability if runs are provided the same seed. */
    srand(seed);
    printf("[!] seed: %ld, num_threads: %d\n", seed, num_threads);
    
    /* Allocate space to populate requests array to be shared by threads. */
    request *requests = calloc(num_requests, sizeof(request));

    /* Initialize requests. */
    for (int i=0; i<num_requests; i++) {
        /* Requests have an equal probability of being a read or a write. */
        requests[i].type = rand() % 2;
        requests[i].duration = (float)(rand() % 10000) / 10000.0;
    }

    /* Allocate threads, as well as allocate and initialize LFDLL, mutex, condition variable, and atomics */
    thread_info *t_info[num_threads];
    pthread_t *threads = (pthread_t*) calloc(num_threads, sizeof(pthread_t));
    
    lfdll_list_t *lf_queue = calloc(1, sizeof(lfdll_list_t));
    LFDLL_INIT(lf_queue, offsetof(queue_element_t, next), offsetof(queue_element_t, prev));
    pthread_mutex_t *lock = calloc(1, sizeof(pthread_mutex_t));
    pthread_mutex_init(lock, NULL);
    pthread_cond_t *cond = calloc(1, sizeof(pthread_cond_t));
    pthread_cond_init(cond, NULL);
    pthread_mutex_t *write_lock = calloc(1, sizeof(pthread_mutex_t));
    pthread_mutex_init(write_lock, NULL);
    pthread_cond_t *write_cond = calloc(1, sizeof(pthread_cond_t));
    pthread_cond_init(write_cond, NULL);
    
    _Atomic int *write = calloc(1, sizeof(_Atomic int));
    atomic_init(write, 0);

    /* Populate the info to be sent to each thread. */
    for (int ts=0; ts<num_threads; ts++) {
        thread_info *info = calloc(1, sizeof(thread_info));
        info->id = ts;
        info->lock = lock;
        info->cond = cond;
        info->write_lock = write_lock;
        info->write_cond = write_cond;
        info->write = write;
        info->queue = lf_queue;
        info->requests = requests;
        info->num_threads = num_threads;
        info->num_requests = num_requests;
        info->verbose = verbose;

        t_info[ts] = info;

        /* Initialize threads with their info */
        pthread_create(&threads[ts], NULL, thread_worker, info);
    }

    /* Collect threads as they terminate. If a thread returns NULL, bail. */
    for (int ts=0; ts<num_threads; ts++) {
        void *ret;
        pthread_join(threads[ts], &ret);

        if (ret == NULL) {
            fprintf(stderr, "[-] Error occured on thread %d\n", ts);
            return -1;
        }
    }

    /* Garbage collection. */
    for (int ts=0; ts<num_threads; ts++) {
        free(t_info[ts]);
    }

    free(requests);
    free(threads);
    free(cond);
    free(lock);
    free(write_lock);
    free(write_cond);
    free(write);

    LFDLL_DESTROY(lf_queue);

    printf("[+] Test success.\n");

    return 0;
}
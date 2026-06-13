Shared Chunk Cache - Serial Dispatch Module

Contents
- lfdll.h, lfdll_errors.h : Lock-Free Doubly Linked List macros/dependencies
- lffl.h, lffl_errors.h   : Lock-Free Free List macros/dependencies

Build test:
  gcc -std=c11 -O2 -Wall -Wextra -o dispatch_test  dispatch_test.c -latomic

Run test:
  ./dispatch_test [-t <num_threads>] [-n <num_requests>] [-s <seed>] [-v]
    - num_threads and num_requests must be greater than 0
    - -v enables verbose output

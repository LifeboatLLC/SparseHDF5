GCC=/usr/bin/gcc
GCC_ARGS=-std=c11 -O2 -Wall -Wextra
L_ARGS=-latomic

run: build
	./dispatch_test -t 32 -n 1000

build:
	${GCC} ${GCC_ARGS} -o disptach_test dispatch_test.c ${L_ARGS}
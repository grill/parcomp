all:
	mpicc -Wall -std=c99 -o stencil-computation stencil-computation.c

debug:
	mpicc -g -Wall -std=c99 -o stencil-computation stencil-computation.c -DDEBUG

optimized:
	mpicc -Wall -std=c99 -o stencil-computation stencil-computation.c -O3

clean:
	rm stencil-computation

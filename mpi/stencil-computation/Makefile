all:
	mpicc -Wall -std=c99 -o stencil-computation stencil-computation.c

debug:
	mpicc -g -Wall -std=c99 -o stencil-computation stencil-computation.c -DDEBUG

debugdebug:
	mpicc -g -Wall -std=c99 -o stencil-computation stencil-computation.c -DDEBUG -DDEBUGDEBUG

optimized:
	mpicc -Wall -std=c99 -o stencil-computation stencil-computation.c -O3 -DEVAL

clean:
	rm stencil-computation

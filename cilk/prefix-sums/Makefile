# Gabriel Grill(1025120)

all: prefixsum

prefixsum: prefixsum.c util.h
	clikc -std=c99 -Wall -fopenmp -lm -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o prefixsum -O3 prefixsum.c

clean:
	rm -rf recursive; rm -rf iterative; rm -rf hillis; rm -rf prefixsum

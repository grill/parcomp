# Gabriel Grill(1025120)

all: prefixsum

prefixsum: prefixsum.c util.h
	gcc -std=c99 -Wall -fopenmp -lm -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o prefixsum -O3 prefixsum.c

iterative: iterative.c util.h
	gcc -std=c99 -Wall -fopenmp -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o iterative -O3 iterative.c

hillis: hillis.c util.h
	gcc -std=c99 -Wall -fopenmp -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o hillis -O3 hillis.c

clean:
	rm -rf recursive; rm -rf iterative; rm -rf hillis; rm -rf prefixsum

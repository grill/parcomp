# Gabriel Grill(1025120)

all: recursive iterative hillis

recursive: recursive.c util.h
	gcc -std=c99 -Wall -fopenmp -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o recursive recursive.c

iterative: iterative.c util.h
	gcc -std=c99 -Wall -fopenmp -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o iterative iterative.c

hillis: hillis.c util.h
	gcc -std=c99 -Wall -fopenmp -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o hillis hillis.c

clean:
	rm -rf recursive; rm -rf iterative; rm -rf hillis;
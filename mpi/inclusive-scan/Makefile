all:
	mpicc -Wall -std=c99 -o inclusive-scan inclusive-scan.c

debug:
	mpicc -g -Wall -std=c99 -o inclusive-scan inclusive-scan.c -DDEBUG

optimized:
	mpicc -Wall -std=c99 -o inclusive-scan inclusive-scan.c -O3 -DEVAL

clean:
	rm inclusive-scan

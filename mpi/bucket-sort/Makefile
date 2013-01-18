# Gabriel Grill(1025120)

all: bucket-sort

run: bucket-sort
	mpirun -np 4 ./bucket-sort

single: bucket-sort
	mpirun -np 1 ./bucket-sort

two: bucket-sort
	mpirun -np 2 ./bucket-sort

radix: bucket-sort
	mpirun -np 4 ./bucket-sort -r

bucket-sort: bucket-sort.c
	mpicc -Wall -lm -o bucket-sort -O3 bucket-sort.c

clean:
	rm -rf bucket-sort

deploy:
	scp * jupiter:/home/ws2012/e1025120

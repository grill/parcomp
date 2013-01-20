# Gabriel Grill(1025120)

all: bucket-sort

run: bucket-sort
	mpirun -np 4 ./bucket-sort

bucket_all: bucket-sort
	mpirun -node 0-35 -nnp 16 bucket-sort -R 100 -n  2147483647

radix_all: bucket-sort
	mpirun -node 0-35 -nnp 16 bucket-sort -r -R 100 -n  2147483647

single: bucket-sort
	mpirun -np 1 ./bucket-sort

two: bucket-sort
	mpirun -np 2 ./bucket-sort

radix: bucket-sort
	mpirun -np 4 ./bucket-sort -r

bucket-sort: bucket-sort.c
	mpicc -Wall -lm -o bucket-sort -O3 bucket-sort.c

radix2: bucket-sort
	mpirun -np 4 ./bucket-sort -R 20 -r -d -n 100

clean:
	rm -rf bucket-sort

deploy:
	scp bucket-sort.c Makefile gendata.sh jupiter:
	ssh jupiter

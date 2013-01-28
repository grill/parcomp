# Gabriel Grill(1025120)

all: prefixsum

prefixsum: prefixsum.c hillis-steele.c recursive.c iterative.c
	gcc -std=c99 -Wall -fopenmp -lm -g -pedantic -DENDEBUG -D_BSD_SOURCE -D_XOPEN_SOURCE=500 -o prefixsum -O3 prefixsum.c hillis-steele.c iterative.c recursive.c

gen_data: recursive_time recursive_work hillis_time hillis_work iterative_time iterative_work reduce

recursive_time:
	./prefixsum -s -p -a -r -n 1000000000 > data/recursive_time_1000000000
	./prefixsum -s -p -a -r -n 500000000 > data/recursive_time_500000000
	./prefixsum -s -p -a -r -n 100000000 > data/recursive_time_100000000
	./prefixsum -s -p -a -r -n 10000000 > data/recursive_time_10000000
	./prefixsum -s -p -a -r -n 1000000 > data/recursive_time_1000000

recursive_work:
	./prefixsum -s -w -a -r -n 100000000 > data/recursive_work_100000000
	./prefixsum -s -w -a -r -n 10000000 > data/recursive_work_10000000
	./prefixsum -s -w -a -r -n 1000000 > data/recursive_work_1000000

hillis_time:
	./prefixsum -s -p -a -h -n 1000000000 > data/hillis_time_1000000000
	./prefixsum -s -p -a -h -n 500000000 > data/hillis_time_500000000
	./prefixsum -s -p -a -h -n 100000000 > data/hillis_time_100000000
	./prefixsum -s -p -a -h -n 10000000 > data/hillis_time_10000000
	./prefixsum -s -p -a -h -n 1000000 > data/hillis_time_1000000

hillis_work:
	./prefixsum -s -w -a -h -n 100000000 > data/hillis_work_100000000
	./prefixsum -s -w -a -h -n 10000000 > data/hillis_work_10000000
	./prefixsum -s -w -a -h -n 1000000 > data/hillis_work_1000000

iterative_time:
	./prefixsum -s -p -a -i -n 1000000000 > data/iterative_time_1000000000
	./prefixsum -s -p -a -i -n 500000000 > data/iterative_time_500000000
	./prefixsum -s -p -a -i -n 100000000 > data/iterative_time_100000000
	./prefixsum -s -p -a -i -n 10000000 > data/iterative_time_10000000
	./prefixsum -s -p -a -i -n 1000000 > data/iterative_time_1000000

iterative_work:
	./prefixsum -s -w -a -i -n 100000000 > data/iterative_work_100000000
	./prefixsum -s -w -a -i -n 10000000 > data/iterative_work_10000000
	./prefixsum -s -w -a -i -n 1000000 > data/iterative_work_1000000

reduce:
	./prefixsum -s -p -a -e -n 1000000000 > data/reduce_1000000000
	./prefixsum -s -p -a -e -n 500000000 > data/reduce_500000000
	./prefixsum -s -p -a -e -n 100000000 > data/reduce_100000000
	./prefixsum -s -p -a -e -n 10000000 > data/reduce_10000000
	./prefixsum -s -p -a -e -n 1000000 > data/reduce_1000000

clean:
	rm -rf prefixsum

deploy:
	scp * saturn:/home/ws2012/e1025120/prefixsum
	ssh saturn
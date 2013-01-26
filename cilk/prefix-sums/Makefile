# Gabriel Grill(1025120)

all: prefixsum

prefixsum: prefixsum.cilk
	/opt/cilk-5.4.6/bin/cilkc -Wall -O3 prefixsum.cilk -o prefixsum

clean:
	rm -rf prefixsum

deploy:
	scp gendata.sh Makefile prefixsum.cilk saturn:/home/ws2012/e1025120/cilk

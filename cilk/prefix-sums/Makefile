# Gabriel Grill(1025120)

all: prefixsum

prefixsum: prefixsum.cilk util.h
	/opt/cilk-5.4.6/bin/cilkc -O3 prefixsum.cilk -o prefixsum

clean:
	rm -rf prefixsum

deploy:
	scp * saturn:/home/ws2012/e1025120/cilk

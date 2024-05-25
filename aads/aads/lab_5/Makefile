all: bench

bench: bench.o matrix.o
	gcc -o bench bench.o matrix.o -lm

matrix.o: matrix.c matrix.h
	gcc -c matrix.c

bench.o: bench.c matrix.h
	gcc -c bench.c

clean:
	rm -rf *.o bench

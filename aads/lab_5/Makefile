W_COMPILE_FLAGS = -Wall -Wextra -Werror

all: bench create_matrix

bench: bench.o matrix.o
	gcc $(W_COMPILE_FLAGS) -o bench bench.o matrix.o -g -lm

matrix.o: matrix.c matrix.h
	gcc $(W_COMPILE_FLAGS) -c matrix.c

bench.o: bench.c matrix.h
	gcc $(W_COMPILE_FLAGS) -c bench.c

create_matrix: create_matrix.c	
	gcc $(W_COMPILE_FLAGS) -o create_matrix create_matrix.c

clean:
	rm -rf *.o bench

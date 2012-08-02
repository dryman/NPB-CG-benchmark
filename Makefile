CC = /usr/local/bin/gcc
CFLAGS = -O3 -Wall -fopenmp
LFLAGS = -Wall -fopenmp
OBJS = cg.o c_print_results.o c_randdp.o c_timers.o wtime.o

all: ${OBJS}
	@echo "Building.."
	${CC} ${CFLAGS} ${OBJS} -o cg

%.o: %.c
	${CC} ${OPTIONS} -c $*.c

clean:
	@echo "Cleaning up *.o"
	-rm -rf *.o
	-rm cg


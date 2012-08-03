CC = /usr/local/bin/gcc
CFLAGS = -g -Wall
LFLAGS = -Wall
OBJS = cg-modified.o c_print_results.o c_randdp.o c_timers.o wtime.o


all: ${OBJS}
	@echo "Building.."
	${CC} ${CFLAGS} ${OBJS} -o cg

run: cg
	./cg

%.o: %.c
	${CC} ${CFLAGS} -c $*.c

clean:
	@echo "Cleaning up *.o"
	-rm -rf *.o
	-rm cg


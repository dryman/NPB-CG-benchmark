CC = /usr/local/bin/gcc
CFLAGS = -g -Wall -fopenmp
LFLAGS = -Wall
OBJS = cg-modified.o c_print_results.o c_randdp.o c_timers.o wtime.o gen-matrix.o

S_FLAGS = -DNA=1400 \
	  -DNONZER=7 \
	  -DNITER=15 \
	  -DSHIFT=10.0 -DRCOND=1.0e-1 -DCONVERTDOUBLE=FALSE \
	  -DCOMPILETIME="03 Aug 2012" \
	  -DNPBVERSION="2.3" \
	  -DCS1=${CC} \
	  -DCS2=${CC} \
	  -DCS3="" -DCS4="" -DCS5=${CFLAGS} -DCS6=${LFLAGS} -DCS7="randdp"

all: ${OBJS}
	@echo "Building.."
	${CC} ${CFLAGS} ${OBJS} -o cg

run: cg
	./cg

%.o: %.c
	${CC} ${CFLAGS} ${S_FLAGS} -c $*.c

clean:
	@echo "Cleaning up *.o"
	-rm -rf *.o
	-rm cg


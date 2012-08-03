#CC = /usr/local/bin/gcc
CC = clang
PARAM_FLAGS = ${CFLAGS} ${A_FLAGS} -DCACHE_LINE_SIZE=${CACHE_LINE_SIZE}
CFLAGS = -O3 -Wall -fopenmp -fblocks
LFLAGS = -Wall
OBJS = c_print_results.o c_randdp.o c_timers.o wtime.o
CFILES = c_print_results.c c_randdp.c c_timers.c wtime.c

CACHE_LINE_SIZE = $(shell sysctl -n hw.cachelinesize)

all: ${OBJS} gen-matrix cg
	@echo "Building.."
	${CC} ${CFLAGS} *.o -o cg

gen-matrix: gen-matrix.c
	${CC} ${PARAM_FLAGS} -c $?

cg: cg-dispatch.c
	${CC} ${PARAM_FLAGS} -c $?

${OBJS}: ${CFILES}
	${CC} ${CFLAGS} -c ${CFILES}

clean:
	@echo "Cleaning up *.o"
	-rm -rf *.o
	-rm cg

S_FLAGS = -DNA=1400 \
	  -DNONZER=7 \
	  -DNITER=15 \
	  -DSHIFT=10.0 -DRCOND=1.0e-1 -DCONVERTDOUBLE=FALSE \
	  -DCOMPILETIME="03 Aug 2012" \
	  -DNPBVERSION="2.3" \
	  -DCS1=${CC} \
	  -DCS2=${CC} \
	  -DCS3="" -DCS4="" -DCS5=${CFLAGS} -DCS6=${LFLAGS} -DCS7="randdp"

W_FLAGS = -DNA=7000 \
	  -DNONZER=8 \
	  -DNITER=15 \
	  -DSHIFT=12.0 -DRCOND=1.0e-1 -DCONVERTDOUBLE=FALSE \
	  -DCOMPILETIME="03 Aug 2012" \
	  -DNPBVERSION="2.3" \
	  -DCS1=${CC} \
	  -DCS2=${CC} \
	  -DCS3="" -DCS4="" -DCS5=${CFLAGS} -DCS6=${LFLAGS} -DCS7="randdp"

A_FLAGS = -DNA=14000 \
	  -DNONZER=11 \
	  -DNITER=15 \
	  -DSHIFT=20.0 -DRCOND=1.0e-1 -DCONVERTDOUBLE=FALSE \
	  -DCOMPILETIME="03 Aug 2012" \
	  -DNPBVERSION="2.3" \
	  -DCS1=${CC} \
	  -DCS2=${CC} \
	  -DCS3="" -DCS4="" -DCS5=${CFLAGS} -DCS6=${LFLAGS} -DCS7="randdp"

B_FLAGS = -DNA=75000 \
	  -DNONZER=13 \
	  -DNITER=75 \
	  -DSHIFT=60.0 -DRCOND=1.0e-1 -DCONVERTDOUBLE=FALSE \
	  -DCOMPILETIME="03 Aug 2012" \
	  -DNPBVERSION="2.3" \
	  -DCS1=${CC} \
	  -DCS2=${CC} \
	  -DCS3="" -DCS4="" -DCS5=${CFLAGS} -DCS6=${LFLAGS} -DCS7="randdp"

C_FLAGS = -DNA=150000 \
	  -DNONZER=15 \
	  -DNITER=75 \
	  -DSHIFT=110.0 -DRCOND=1.0e-1 -DCONVERTDOUBLE=FALSE \
	  -DCOMPILETIME="03 Aug 2012" \
	  -DNPBVERSION="2.3" \
	  -DCS1=${CC} \
	  -DCS2=${CC} \
	  -DCS3="" -DCS4="" -DCS5=${CFLAGS} -DCS6=${LFLAGS} -DCS7="randdp"


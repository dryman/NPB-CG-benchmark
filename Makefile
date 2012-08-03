CC = /usr/local/bin/gcc
CFLAGS = -O3 -Wall -fopenmp
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

all: ${OBJS}
	@echo "Building.."
	${CC} ${CFLAGS} ${OBJS} -o cg

run: cg
	./cg

%.o: %.c
	${CC} ${CFLAGS} ${B_FLAGS} -c $*.c

clean:
	@echo "Cleaning up *.o"
	-rm -rf *.o
	-rm cg


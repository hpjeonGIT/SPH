#
# Makefile for SPH program SPH2D
# 2001.08.31
.SUFFIXES: .o .f90
#
# Compiler
F90 = f90
#
# Compiling option
FLAGS  = +O2
LIB = -lm
#
# Object files
OBJTS = artvisco.o kernel.o mainsph.o parsor.o solver.o strength.o eos.o
#
# Object file compiling
.f90.o:
	${F90} -c ${FLAGS} $*.f90 ${LIB}
#
# Execution file
SPHax:${OBJTS}
	${F90} -o SPHax ${OBJTS} 
#
# clean
clean:
	rm -rf *.o core

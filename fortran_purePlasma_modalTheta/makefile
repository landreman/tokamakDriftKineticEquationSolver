# Makefile for the local tokamak drift-kinetic equation solver

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

TARGET = tokamakDriftKineticEquationSolver

OBJS = $(TARGET).o makeXGrid.o geometry.o writeOutput.o parameters.o convergenceScan.o readInput.o fft.o \
	modalMultiplicationMatrix.o poldif.o uniformDiffMatrices.o interpolationMatrix.o \
	polyInterp.o solveDKE.o fftpack/fftpack.a quadpack/quadpack.a

all: $(OBJS)
	$(FLINKER) -o $(TARGET) $(OBJS) $(PETSC_LIB)

.PHONY: all cleaner

%.o: %.F90
	$(FLINKER) -c $<

fftpack/fftpack.a:
	cd fftpack; make all

quadpack/quadpack.a:
	cd quadpack; make all

cleaner:
	rm -f *.o *.mod *.a *.dat *~ $(TARGET)
	cd fftpack; rm -f *.o *.a
	cd quadpack; rm -f *.o *.a

# Dependencies:
geometry.o: parameters.o modalMultiplicationMatrix.o
writeOuput.o: parameters.o
readInput.o: parameters.o
modalMultiplicationMatrix.o: fft.o parameters.o
solveDKE.o: parameters.o makeXGrid.o geometry.o modalMultiplicationMatrix.o poldif.o uniformDiffMatrices.o interpolationMatrix.o polyInterp.o
convergenceScan.o: parameters.o
$(TARGET).o: parameters.o geometry.o writeOutput.o solveDKE.o readInput.o convergenceScan.o


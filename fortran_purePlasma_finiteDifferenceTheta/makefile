# Makefile for the local tokamak drift-kinetic solver

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

TARGET = tokamakDriftKineticEquationSolver

OBJS = $(TARGET).o makeXGrid.o geometry.o writeOutput.o parameters.o convergenceScan.o readInput.o \
	poldif.o uniformDiffMatrices.o interpolationMatrix.o \
	polyInterp.o solveDKE.o quadpack/quadpack.a

all: $(OBJS)
	$(FLINKER) -o $(TARGET) $(OBJS) $(PETSC_LIB)

.PHONY: all cleaner

%.o: %.F90
	$(FLINKER) -c $<

quadpack/quadpack.a:
	cd quadpack; make all

cleaner: clean
	rm -f *.o *.mod *.a *.dat *~ $(TARGET)
	cd quadpack; rm -f *.o *.a

# Dependencies:
geometry.o: parameters.o
writeOutput.o: parameters.o convergenceScan.o
readInput.o: parameters.o
solveDKE.o: parameters.o makeXGrid.o geometry.o poldif.o uniformDiffMatrices.o interpolationMatrix.o polyInterp.o
convergenceScan.o: parameters.o
$(TARGET).o: parameters.o geometry.o writeOutput.o solveDKE.o readInput.o convergenceScan.o

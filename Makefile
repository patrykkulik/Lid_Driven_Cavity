##############################################################################
# FILE: Makefile.MPI.c
# DESCRIPTION:
# Makefile for the HPC Lid Driven Cavity code
# AUTHOR: Patryk Kulik
###############################################################################

#MPI compiler:
CXX    =   mpicxx


#Add flags
# CXXFLAGS   =   -O2 -w

#Add libraries to link
LDLIBS = -lscalapack-openmpi -llapack -lblas

all:    LidDrivenCavitySolver

clean:
	/bin/rm -rf     \
	LidDrivenCavitySolver        \
	*.o


LidDrivenCavitySolver:  LidDrivenCavitySolver.cpp
	$(CXX) LidDrivenCavitySolver.cpp $(LDLIBS) -o LidDrivenCavitySolver



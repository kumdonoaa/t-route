# compiler
FC := $(F90)

# compile flags
#FCFLAGS = -c -fdefault-real-8 -fno-align-commons -fbounds-check --free-form
FCFLAGS = -g -c -O2 -fPIC
# link flags
FLFLAGS = -static-gfortran -static-libgcc -no-defaultlibs -lgfortran -lgcc
VPATH = ../Reservoir_singleTS

.PHONY: reach reservoir all install
all: reach reservoir

reach: mc_single_seg.o pymc_single_seg.o

mc_single_seg.o: MCsingleSegStime_f2py_NOLOOP.f90 precis.mod
	$(FC) $(FCFLAGS) -o $@ $<

pymc_single_seg.o: pyMCsingleSegStime_NoLoop.f90
	$(FC) $(FCFLAGS) -o $@ $<

reservoir:  precis.mod module_levelpool.o pymodule_levelpool.o

precis.mod: varPrecision.f90
	$(FC) $(FCFLAGS) $<

install:
	cp *.o ../../../../src/python_routing_v02/troute/routing/fast_reach

%.o: %.F
	$(FC) $(FCFLAGS) -o $@ $<

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod

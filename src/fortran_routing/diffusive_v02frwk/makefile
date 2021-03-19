# compiler
FC := gfortran

# compile flags
FCFLAGS = -g -c -O2 -fPIC
# link flags
# FLFLAGS = -static-gfortran -static-libgcc -no-defaultlibs -lgfortran -lgcc
# VPATH = ../Reservoir_singleTS

#.PHONY: reach reservoir all install
#all: reach reservoir

#reach: mc_single_seg.o pymc_single_seg.o

diffusive.o: maindiff_v2.f90
	$(FC) $(FCFLAGS) -o $@ $<

pydiffusive.o: pymaindiff_v2.f90
	$(FC) $(FCFLAGS) -o $@ $<

install:
	cp *.o ../../../src/python_routing_v02/fast_reach

clean:
	rm -f *.o *.mod

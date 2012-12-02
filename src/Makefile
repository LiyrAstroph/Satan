FC=gfortran -O2 -Wall -m64

all:diskmain speccal

diskmain: diskmain.f90 const.mod diskvars.mod disksol.o diskset.o diskrad.o diskcor.o specfun.o \
	newton.o Makefile
	$(FC) -o $@ $< disksol.o diskset.o diskrad.o diskcor.o specfun.o newton.o

const.mod: const.f90
	$(FC) -fsyntax-only $< -o $@

diskvars.mod: diskvars.f90
	$(FC) -fsyntax-only $< -o $@

disksol.o: disksol.f90 const.mod diskvars.mod
	$(FC) -c $< -o $@

diskset.o: diskset.f90 const.mod diskvars.mod
	$(FC) -c $< -o $@

diskrad.o: diskrad.f90 const.mod diskvars.mod
	$(FC) -c $< -o $@

diskcor.o: diskcor.f90 const.mod diskvars.mod
	$(FC) -c $< -o $@

newton.o: newton.f90 const.mod diskvars.mod
	$(FC) -c $< -o $@

specfun.o: specfun.f90
	$(FC) -c $< -o $@

speccal: speccal.f90 speccomp.o specvars.mod specfun.o phindex.o
	$(FC) -o $@ $< speccomp.o specfun.o phindex.o

specvars.mod: specvars.f90
	$(FC) -fsyntax-only $< -o $@

speccomp.o: speccomp.f90 specvars.mod
	$(FC) -c $< -o $@

phindex.o: phindex.f90 specvars.mod const.mod
	$(FC) -c $< -o $@

clean:
	rm diskmain *.o *.mod

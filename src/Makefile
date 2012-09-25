FC=gfortran -O2

all:diskmain

diskmain: diskmain.f90 const.mod diskvars.mod disksol.o diskset.o Makefile
	$(FC) -o $@ $< disksol.o diskset.o

const.mod: const.f90
	$(FC) -fsyntax-only $< -o $@

diskvars.mod: diskvars.f90
	$(FC) -fsyntax-only $< -o $@

disksol.o: disksol.f90 const.mod diskvars.mod
	$(FC) -c $< -o $@

diskset.o: diskset.f90 const.mod diskvars.mod
	$(FC) -c $< -o $@


clean:
	rm diskmain *.o *.mod

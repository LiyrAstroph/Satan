FC=gfortran -O2

all:diskmain

diskmain: diskmain.f90 const.o diskvars.o disksol.o diskset.o
	$(FC) -o $@ $+

const.o: const.f90
	$(FC) -c $<

diskvars.o: diskvars.f90
	$(FC) -c $<

disksol.o: disksol.f90
	$(FC) -c $<

diskset.o: diskset.f90
	$(FC) -c $<


clean:
	rm diskmain *.o

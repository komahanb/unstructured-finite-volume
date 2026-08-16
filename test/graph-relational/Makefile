include ../../Makefile.in

#=====================================================================#
# The relational view over the production kernel.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: run lifetime refusal

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

lifetime: lifetime.o
	$(F90) lifetime.o -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) refusal.o -o $@ $(LDFLAGS)

print-f90:
	@echo $(F90)

print-std:
	@echo $(FSTD)

clean:
	@rm -f *.o *.mod *.so run lifetime refusal
	@rm -f fortran-assignment/*.o fortran-assignment/*.mod \
	       fortran-assignment/candidate fortran-assignment/probe.out

.SUFFIXES: .f90 .o
.PHONY: all clean print-f90 print-std

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

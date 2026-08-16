include ../../Makefile.in

#=====================================================================#
# Two programs: the law suite, and the refusals that must die.
# -Wtarget-lifetime keeps the borrows honest at compile time.
#=====================================================================#

# Ufvm (-I../../src for scalar.fpp, -I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a -llapack

all: run refusal

run: listed_set_fixture.o test.o
	$(F90) relational_fixture.o listed_set_fixture.o test.o -o $@ $(LDFLAGS)

refusal: relational_fixture.o refusal.o
	$(F90) relational_fixture.o refusal.o -o $@ $(LDFLAGS)

test.o: relational_fixture.o listed_set_fixture.o

listed_set_fixture.o: ../graph-relation/listed_set_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@

refusal.o: relational_fixture.o

relational_fixture.o: ../support/relational_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

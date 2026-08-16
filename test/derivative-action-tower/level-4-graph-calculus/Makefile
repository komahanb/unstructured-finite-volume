include ../../../Makefile.in

#=====================================================================#
# Derivative action tower, Gate A: one rung, its own truths. The
# framework modules this level may see are named in the tower's
# import gate (check_imports.sh) - the gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: derivative_assert.o test.o
	$(F90) relational_fixture.o derivative_assert.o test.o -o $@ $(LDFLAGS)

test.o: relational_fixture.o derivative_assert.o

derivative_assert.o: ../common/derivative_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

relational_fixture.o: ../../support/relational_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

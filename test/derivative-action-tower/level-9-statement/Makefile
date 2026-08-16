include ../../../Makefile.in

#=====================================================================#
# Derivative level 9: the statement suite. The composition rung:
# carriers, relations, algebra, the relational graph, the profile
# and algorithms for the derived order, the field - and the REUSED
# level-8 derivative constitution, compiled from its own rung's
# directory. There is deliberately no production statement module
# to import, no minimization anywhere on this road, and no adapter:
# nothing here must satisfy a legacy operation face.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: derivative_assert.o derivative_constitution_fixture.o test.o
	$(F90) relational_fixture.o derivative_assert.o derivative_constitution_fixture.o test.o -o $@ $(LDFLAGS)

test.o: relational_fixture.o derivative_constitution_fixture.o
derivative_constitution_fixture.o: ../level-8-derivative-constitution/derivative_constitution_fixture.f90 derivative_assert.o
	$(F90) $(F90FLAGS) -c $< -o $@

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

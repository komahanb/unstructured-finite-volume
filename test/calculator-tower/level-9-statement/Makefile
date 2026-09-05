include ../../../Makefile.in

#=====================================================================#
# Level 9: the statement suite. The composition rung: carriers,
# relations, algebra, the relational graph, the field, the legacy
# operation host, minimization with GMRES - and the REUSED level-8
# constitution fixture, compiled from its own rung's directory.
# There is deliberately no production statement module to import.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run

run: calculator_assert.o constitution_fixture.o test.o
	$(F90) calculator_assert.o constitution_fixture.o test.o -o $@ $(LDFLAGS)

test.o: constitution_fixture.o
constitution_fixture.o: ../level-8-constitution/constitution_fixture.f90 calculator_assert.o
	$(F90) $(F90FLAGS) -c $< -o $@

calculator_assert.o: ../common/calculator_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

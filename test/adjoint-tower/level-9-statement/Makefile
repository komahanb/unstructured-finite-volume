include ../../../Makefile.in

#=====================================================================#
# Adjoint level 9: the statement suite. The composition rung: the
# relation stack that carries the structural supports, the model
# graph that owns them, the field, the legacy operation host and
# the GMRES citizen - plus the REUSED level-8 constitution,
# compiled from its own rung's directory. No new law, no new type,
# and deliberately no class_graph_linearization: the statement does
# not retreat into the specialized same-domain path.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run

run: adjoint_assert.o adjoint_constitution_fixture.o test.o
	$(F90) adjoint_assert.o adjoint_constitution_fixture.o test.o -o $@ $(LDFLAGS)

test.o: adjoint_constitution_fixture.o
adjoint_constitution_fixture.o: ../level-8-constitution/adjoint_constitution_fixture.f90 adjoint_assert.o
	$(F90) $(F90FLAGS) -c $< -o $@

adjoint_assert.o: ../common/adjoint_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

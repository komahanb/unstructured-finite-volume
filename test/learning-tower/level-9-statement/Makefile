include ../../../Makefile.in

#=====================================================================#
# Learning level 9: the statement suite. The composition rung:
# sets, relations, algebra, the related graph, the interpretation
# and algorithms for the derived order, the field, the legacy
# operation host, the GMRES citizen - and the REUSED level-8
# constitution fixture, compiled from its own rung's directory.
# There is deliberately no production statement module to import,
# and Level 7's affine oracle appears nowhere in this path.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: learning_assert.o learning_constitution_fixture.o constituted_residual_fixture.o test.o
	$(F90) learning_assert.o learning_constitution_fixture.o constituted_residual_fixture.o test.o -o $@ $(LDFLAGS)

test.o: constituted_residual_fixture.o
constituted_residual_fixture.o: learning_constitution_fixture.o
learning_constitution_fixture.o: ../level-8-constitution/learning_constitution_fixture.f90 learning_assert.o
	$(F90) $(F90FLAGS) -c $< -o $@

learning_assert.o: ../common/learning_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

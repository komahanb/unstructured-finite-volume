include ../../../Makefile.in

#=====================================================================#
# Adjoint level 8: the one-constitution suite. A single coefficient
# table, keyed by member, generates the primal residual, the
# response, and the forward and reverse actions - all driven by the
# Gate-A structural supports. The whole battery runs twice under
# independently permuted role enumerations. The import gate holds
# the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run

run: adjoint_assert.o adjoint_constitution_fixture.o test.o
	$(F90) adjoint_assert.o adjoint_constitution_fixture.o test.o -o $@ $(LDFLAGS)

test.o: adjoint_constitution_fixture.o
adjoint_constitution_fixture.o: adjoint_assert.o

adjoint_assert.o: ../common/adjoint_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

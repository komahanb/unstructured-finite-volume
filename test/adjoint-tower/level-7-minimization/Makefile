include ../../../Makefile.in

#=====================================================================#
# Adjoint level 7: the solver-neutrality suite and the same-size
# refusals. The equations are SUPPLIED here, not generated - Level 8
# removes that independence. gmres inherits attach/constant/apply
# from the minimizer base, so graph_minimization is not imported
# directly. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run refusal

run: adjoint_assert.o opaque_equation_fixture.o test.o
	$(F90) adjoint_assert.o opaque_equation_fixture.o test.o -o $@ $(LDFLAGS)

refusal: adjoint_assert.o opaque_equation_fixture.o refusal.o
	$(F90) adjoint_assert.o opaque_equation_fixture.o refusal.o -o $@ $(LDFLAGS)

test.o: opaque_equation_fixture.o
refusal.o: opaque_equation_fixture.o
opaque_equation_fixture.o: adjoint_assert.o

adjoint_assert.o: ../common/adjoint_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

include ../../../Makefile.in

#=====================================================================#
# Derivative level 8: the numerical action suite and the refusals
# that must die. The constitution is TEST-LOCAL - primal laws and
# one local linearization per operation, port-relative, slot-blind.
# Tangent and reverse are two actions of that one shadow. No
# solver, no minimization, no Jacobian assembly, no binary storage:
# the import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run refusal

run: derivative_assert.o derivative_constitution_fixture.o test.o
	$(F90) derivative_assert.o derivative_constitution_fixture.o test.o -o $@ $(LDFLAGS)

refusal: derivative_assert.o derivative_constitution_fixture.o refusal.o
	$(F90) derivative_assert.o derivative_constitution_fixture.o refusal.o -o $@ $(LDFLAGS)

test.o: derivative_assert.o derivative_constitution_fixture.o
refusal.o: derivative_assert.o derivative_constitution_fixture.o
derivative_constitution_fixture.o: derivative_assert.o

derivative_assert.o: ../common/derivative_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

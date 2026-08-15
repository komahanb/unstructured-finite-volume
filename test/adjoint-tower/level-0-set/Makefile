include ../../../Makefile.in

#=====================================================================#
# Adjoint tower: one rung, its own truths. The framework modules
# this level may see are named in the tower's import gate
# (check_imports.sh) - the gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: adjoint_assert.o test.o
	$(F90) adjoint_assert.o test.o -o $@ $(LDFLAGS)

test.o: adjoint_assert.o

adjoint_assert.o: ../common/adjoint_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

include ../../../Makefile.in

#=====================================================================#
# Learning level 7: the fitting suite. Framework modules this level
# may see: the carrier/domain machinery, the field, the legacy
# graph_operation host contract and the GMRES citizen - which
# inherits attach, constant and solve from the minimizer base, so
# graph_minimization itself stays off the allowlist. Constitution
# and statement remain forbidden: the solver never learns what
# predict or error mean. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: learning_assert.o test.o
	$(F90) learning_assert.o test.o -o $@ $(LDFLAGS)

test.o: learning_assert.o

learning_assert.o: ../common/learning_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

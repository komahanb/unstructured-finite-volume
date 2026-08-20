include ../../../Makefile.in

#=====================================================================#
# Level 7: the solve suite. Framework modules this level may see:
# the carrier/domain machinery, the field, the legacy
# graph_operation host contract, graph_minimization and the GMRES
# citizen. Constitution and statement remain forbidden: the solver
# never learns what + or x mean. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run

run: calculator_assert.o test.o
	$(F90) calculator_assert.o test.o -o $@ $(LDFLAGS)

test.o: calculator_assert.o

calculator_assert.o: ../common/calculator_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

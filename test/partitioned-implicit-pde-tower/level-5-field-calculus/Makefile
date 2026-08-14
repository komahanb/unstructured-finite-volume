include ../../../Makefile.in

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: partitioned_pde_assert.o test.o
	$(F90) partitioned_pde_assert.o test.o -o $@ $(LDFLAGS)

test.o: partitioned_pde_assert.o
partitioned_pde_assert.o: ../common/partitioned_pde_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

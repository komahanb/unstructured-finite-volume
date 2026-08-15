include ../../../Makefile.in

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: time_assert.o time_sets_fixture.o test.o
	$(F90) time_assert.o time_sets_fixture.o test.o -o $@ $(LDFLAGS)

test.o: time_assert.o time_sets_fixture.o
time_sets_fixture.o: ../common/time_sets_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@
time_assert.o: ../common/time_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

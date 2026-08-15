include ../../../Makefile.in

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: visualization_assert.o visualization_sets_fixture.o test.o
	$(F90) visualization_assert.o visualization_sets_fixture.o test.o -o $@ $(LDFLAGS)

test.o: visualization_assert.o visualization_sets_fixture.o
visualization_sets_fixture.o: ../common/visualization_sets_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@
visualization_assert.o: ../common/visualization_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

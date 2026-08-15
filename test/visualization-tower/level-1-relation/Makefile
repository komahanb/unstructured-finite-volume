include ../../../Makefile.in

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

FIXTURES = visualization_assert.o visualization_sets_fixture.o \
           visualization_relations_fixture.o

all: run

run: $(FIXTURES) test.o
	$(F90) $(FIXTURES) test.o -o $@ $(LDFLAGS)

test.o: $(FIXTURES)
visualization_relations_fixture.o: ../common/visualization_relations_fixture.f90 visualization_assert.o
	$(F90) $(F90FLAGS) -c $< -o $@
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

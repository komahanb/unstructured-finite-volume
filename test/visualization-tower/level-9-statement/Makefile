include ../../../Makefile.in

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

FIXTURES = visualization_assert.o visualization_carriers_fixture.o \
           visualization_relations_fixture.o visualization_algebra_fixture.o \
           structural_renderer_fixture.o visualization_values_fixture.o \
           valued_renderer_fixture.o production_pattern_renderer_fixture.o

all: run

run: $(FIXTURES) test.o
	$(F90) $(FIXTURES) test.o -o $@ $(LDFLAGS)

test.o: $(FIXTURES)
production_pattern_renderer_fixture.o: ../common/production_pattern_renderer_fixture.f90 structural_renderer_fixture.o visualization_carriers_fixture.o
	$(F90) $(F90FLAGS) -c $< -o $@
valued_renderer_fixture.o: ../common/valued_renderer_fixture.f90 structural_renderer_fixture.o visualization_carriers_fixture.o
	$(F90) $(F90FLAGS) -c $< -o $@
visualization_values_fixture.o: ../common/visualization_values_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@
structural_renderer_fixture.o: ../common/structural_renderer_fixture.f90 visualization_carriers_fixture.o
	$(F90) $(F90FLAGS) -c $< -o $@
visualization_algebra_fixture.o: ../common/visualization_algebra_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@
visualization_relations_fixture.o: ../common/visualization_relations_fixture.f90 visualization_assert.o
	$(F90) $(F90FLAGS) -c $< -o $@
visualization_carriers_fixture.o: ../common/visualization_carriers_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@
visualization_assert.o: ../common/visualization_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

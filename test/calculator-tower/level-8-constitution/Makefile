include ../../../Makefile.in

#=====================================================================#
# Level 8: the meaning suite. Framework modules this level may see:
# the carriers/domains, the relations, the algebra, the field, and
# the TEST-LOCAL arithmetic constitution fixture beside this file.
# No solver, no graph host, no graph grammar, no statement: the
# import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run

run: calculator_assert.o constitution_fixture.o test.o
	$(F90) calculator_assert.o constitution_fixture.o test.o -o $@ $(LDFLAGS)

test.o: constitution_fixture.o
constitution_fixture.o: calculator_assert.o

test.o: calculator_assert.o

calculator_assert.o: ../common/calculator_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

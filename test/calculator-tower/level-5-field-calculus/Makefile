include ../../../Makefile.in

#=====================================================================#
# Level 5: the field suite. The view and algorithm refusal
# laws are pinned generically in test/graph-algorithms; this level
# proves the walk. Framework modules this level may see: the
# carriers, the relations, the algebra, the container, the profile
# needs a domain, not a graph - the import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

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

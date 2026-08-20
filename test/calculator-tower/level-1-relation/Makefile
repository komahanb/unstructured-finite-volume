include ../../../Makefile.in

#=====================================================================#
# Level 1: the law suite and the refusals that must die. The
# framework modules this level may see: the set modules and
# graph_relation - the import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run refusal

run: calculator_assert.o test.o
	$(F90) calculator_assert.o test.o -o $@ $(LDFLAGS)

refusal: calculator_assert.o refusal.o
	$(F90) calculator_assert.o refusal.o -o $@ $(LDFLAGS)

test.o: calculator_assert.o
refusal.o: calculator_assert.o

calculator_assert.o: ../common/calculator_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

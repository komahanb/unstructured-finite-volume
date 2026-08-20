include ../../../Makefile.in

#=====================================================================#
# Learning level 1: the law suite and the refusals that must die.
# The framework modules this level may see: the set modules and
# relation_finitary - no graph concept compiles here. The import gate
# holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run refusal

run: learning_assert.o test.o
	$(F90) learning_assert.o test.o -o $@ $(LDFLAGS)

refusal: learning_assert.o refusal.o
	$(F90) learning_assert.o refusal.o -o $@ $(LDFLAGS)

test.o: learning_assert.o
refusal.o: learning_assert.o

learning_assert.o: ../common/learning_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

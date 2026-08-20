include ../../../Makefile.in

#=====================================================================#
# Learning level 5: the value suite. A field needs a domain, not a
# graph: this rung sees learning_assert, the set modules and
# class_graph_field - no relation, no container, no profile, no
# minimizer. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run

run: learning_assert.o test.o
	$(F90) learning_assert.o test.o -o $@ $(LDFLAGS)

test.o: learning_assert.o

learning_assert.o: ../common/learning_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

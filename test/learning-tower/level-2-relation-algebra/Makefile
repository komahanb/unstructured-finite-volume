include ../../../Makefile.in

#=====================================================================#
# Learning level 2: the derivation suite. The generic algebra laws
# live in test/graph-algebra; this level proves THIS structure
# derives THIS dependency. Framework modules this level may see:
# graph_set, graph_relation, graph_relation_algebra - D is held
# as class(relation), so the storage behind the composition is
# never named. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

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

include ../../../Makefile.in

#=====================================================================#
# Level 6: the sparsity suite. Framework modules this level may
# see: graph_carrier, graph_relation, graph_relation_algebra, and
# graph_binary_relation - the last DELIBERATELY, because the
# reverse structure is transpose_of over the derived J, part of the
# level's dependency boundary. The import gate holds the door.
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

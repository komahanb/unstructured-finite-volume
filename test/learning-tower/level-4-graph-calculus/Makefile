include ../../../Makefile.in

#=====================================================================#
# Learning level 4: the interpretation suite. The view and
# algorithm refusal laws are pinned generically in
# test/graph-algorithms; this rung proves THIS derived dependency
# can be read and walked. Framework modules this level may see: the
# carriers, the relations, the algebra, the container, the profile
# and the algorithms - never the binary storage the view leans on
# internally. The import gate holds the door.
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

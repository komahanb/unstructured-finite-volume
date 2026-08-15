include ../../../Makefile.in

#=====================================================================#
# Learning level 3: the container suite and the refusals that must
# die. The main test sees the sets, the relations, the algebra
# and the container; graph_binary_relation is granted for the VIEW
# refusal alone (the production path never touches it). The import
# gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

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

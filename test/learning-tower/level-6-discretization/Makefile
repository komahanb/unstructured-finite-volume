include ../../../Makefile.in

#=====================================================================#
# Learning level 6: the structural-Jacobian rung. The binary
# citizen is earned at last - J_Theta must materialize its
# generated pairs and answer its reverse as a transpose view of the
# same stored truth. Fields stay forbidden: dependency structure
# belongs to the model, not to one data instance. The import gate
# holds the door.
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

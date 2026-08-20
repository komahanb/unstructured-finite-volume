include ../../../Makefile.in

#=====================================================================#
# Adjoint level 4: the implicit-coupling suite and the refusal that
# must die. What is proved here is negative - an implicitly coupled
# state system has no topological order - so the refusal binary
# carries the rung's central truth. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a

all: run refusal

run: adjoint_assert.o test.o
	$(F90) adjoint_assert.o test.o -o $@ $(LDFLAGS)

refusal: adjoint_assert.o refusal.o
	$(F90) adjoint_assert.o refusal.o -o $@ $(LDFLAGS)

test.o: adjoint_assert.o
refusal.o: adjoint_assert.o

adjoint_assert.o: ../common/adjoint_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

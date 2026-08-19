include ../../Makefile.in

#=====================================================================#
# Two programs over the dense direct seat: the law suite proving
# one more concretion of the graph minimization tower - matrices
# as stencil graphs, elimination behind the minimizer face - and
# the refusals that must die. -Wtarget-lifetime keeps the borrows
# honest at compile time.
#=====================================================================#

# Ufvm (-I../../src for scalar.fpp, -I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a -llapack

all: run refusal

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) refusal.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

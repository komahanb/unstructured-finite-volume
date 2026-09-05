include ../../Makefile.in

#=====================================================================#
# Two-type graph kernel, the views over it, the test suite, and the
# refusals. Links libufvm only for token_identity's token.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: run refusal mutation

# graph_fractal now lives in src and arrives through libufvm.a; this
# directory builds only the views and the programs over them.

run: graph_views.o test.o
	$(F90) graph_views.o test.o -o $@ $(LDFLAGS)

refusal: graph_views.o refusal.o
	$(F90) graph_views.o refusal.o -o $@ $(LDFLAGS)

mutation: graph_views.o mutation.o
	$(F90) graph_views.o mutation.o -o $@ $(LDFLAGS)

test.o: graph_views.o
refusal.o: graph_views.o
mutation.o: graph_views.o

clean:
	@rm -f *.o *.mod *.so run refusal mutation
	@rm -f fortran-recursion/*.mod fortran-mutation/*.mod

# run.sh compiles the declaration-order probes with this compiler and
# this standard.
print-f90:
	@echo $(F90)

print-std:
	@echo $(FSTD)

.SUFFIXES: .f90 .o
.PHONY: all clean print-f90 print-std

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

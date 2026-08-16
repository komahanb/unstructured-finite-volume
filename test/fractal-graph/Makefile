include ../../Makefile.in

#=====================================================================#
# Two-type graph kernel, the views over it, the test suite, and the
# refusals. Links libufvm only for graph_identity's token.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: run refusal

run: fractal_graph.o graph_views.o test.o
	$(F90) fractal_graph.o graph_views.o test.o -o $@ $(LDFLAGS)

refusal: fractal_graph.o graph_views.o refusal.o
	$(F90) fractal_graph.o graph_views.o refusal.o -o $@ $(LDFLAGS)

test.o: graph_views.o fractal_graph.o
refusal.o: graph_views.o fractal_graph.o
graph_views.o: fractal_graph.o

clean:
	@rm -f *.o *.mod *.so run refusal fortran-recursion/*.mod

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

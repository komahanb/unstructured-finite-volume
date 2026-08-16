include ../../Makefile.in

#=====================================================================#
# The spike: kernel, views, the law suite, and the refusals that
# must die. Links libufvm only for graph_identity's token.
#=====================================================================#

# Ufvm (-I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a -llapack

all: run refusal

run: fractal_graph.o fractal_views.o test.o
	$(F90) fractal_graph.o fractal_views.o test.o -o $@ $(LDFLAGS)

refusal: fractal_graph.o fractal_views.o refusal.o
	$(F90) fractal_graph.o fractal_views.o refusal.o -o $@ $(LDFLAGS)

test.o: fractal_views.o fractal_graph.o
refusal.o: fractal_views.o fractal_graph.o
fractal_views.o: fractal_graph.o

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

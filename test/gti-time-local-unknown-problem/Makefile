include ../../Makefile.in

#=====================================================================#
# Two programs over the unknown seat: the law suite proving
# injection and residual(q_trial), and the refusals that must die
# - including a mock foreign field the mint must refuse.
# -Wtarget-lifetime keeps the borrows honest at compile time.
#=====================================================================#

# Ufvm (-I../../src for scalar.fpp, -I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a -llapack

all: run refusal

run: toy_forms.o test.o
	$(F90) toy_forms.o test.o -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) refusal.o -o $@ $(LDFLAGS)

test.o: toy_forms.o

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

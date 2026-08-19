include ../../Makefile.in

#=====================================================================#
# Two programs over the controller seat: the law suite proving the
# propose/try/measure/judge/retry loop against closed BDF roots -
# every decision the policy's, every execution the mechanisms',
# not one driver line changed - and the refusals that must die
# when a policy contradicts itself. -Wtarget-lifetime keeps the
# borrows honest at compile time.
#=====================================================================#

# Ufvm (-I../../src for scalar.fpp, -I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a -llapack

all: run refusal

run: toy_forms.o test.o
	$(F90) toy_forms.o test.o -o $@ $(LDFLAGS)

refusal: toy_forms.o refusal.o
	$(F90) toy_forms.o refusal.o -o $@ $(LDFLAGS)

test.o: toy_forms.o
refusal.o: toy_forms.o

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

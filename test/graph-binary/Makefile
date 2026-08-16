include ../../Makefile.in

#=====================================================================#
# Two programs: the law suite, and the refusals that must die. The
# listed-set fixture is borrowed from the relation suite - one
# fixture, one truth.
#=====================================================================#

# Ufvm (-I../../src for scalar.fpp, -I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/
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

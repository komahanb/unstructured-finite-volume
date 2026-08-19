include ../../Makefile.in

#=====================================================================#
# Two programs over the attached-value change seat: the law suite
# proving one concrete value change under the ONE generic
# controller - restoring UNATTACHED, UNKNOWN, and KNOWN prior
# states exactly - and the refusals that must die at their owning
# seats. -Wtarget-lifetime keeps the borrows honest at compile
# time.
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

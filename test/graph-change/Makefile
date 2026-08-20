include ../../Makefile.in

#=====================================================================#
# Builds two binaries: run (the assertion suite) and refusal (the
# invalid-input cases). Both link the shared fixture
# toy_changes.o. -Wtarget-lifetime warns when a pointer may
# outlive its target.
#=====================================================================#

# Ufvm (-I../../src for scalar.fpp, -I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a

all: run refusal

run: toy_changes.o test.o
	$(F90) toy_changes.o test.o -o $@ $(LDFLAGS)

refusal: toy_changes.o refusal.o
	$(F90) toy_changes.o refusal.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

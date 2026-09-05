include ../../Makefile.in

#=====================================================================#
# Builds two binaries: run (the assertion suite) and refusal (the
# invalid-input cases). Both link the shared fixture toy_forms.o.
# -Wtarget-lifetime warns when a pointer may outlive its target.
#=====================================================================#

# Ufvm (-I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a

all: run refusal

run: toy_forms.o test.o
	$(F90) toy_forms.o test.o -o $@ $(LDFLAGS)

refusal: toy_forms.o refusal.o
	$(F90) toy_forms.o refusal.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

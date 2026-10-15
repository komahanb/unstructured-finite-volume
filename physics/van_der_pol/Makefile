include ../../Makefile.in

#=====================================================================#
# THE VAN DER POL OSCILLATOR, stated once: the physics, the adjoint,
# the sensitivity along the design and the design. src/ is the
# mathematics and the computing; the program here states the fields,
# the laws with their data, the objective and the Lagrangian, and
# calls minimize.
# -Wtarget-lifetime checks pointer target lifetimes at compile time.
#=====================================================================#

# Ufvm (-I../../lib for the .mod files; lib_par for COARRAY=lib,
# whose link line Makefile.in supplies)
LIB = ../../lib
VARIANT = serial
ifeq ($(COARRAY),lib)
  LIB = ../../lib_par
  VARIANT = coarray
endif
F90FLAGS += -I$(LIB)/ -Wtarget-lifetime
LDFLAGS += -L$(LIB)/ $(LIB)/libufvm.a

PROGRAMS = van_der_pol

all: $(PROGRAMS)

$(PROGRAMS): %: %.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so $(PROGRAMS) variant.stamp

.SUFFIXES: .f90 .o
.PHONY: all clean FORCE

# an object depends on the library too: a rebuilt library may have
# changed a module interface the program compiled against; and on
# the build variant, so that a serial build after a coarray one (or
# the reverse) recompiles rather than relinks the other's objects
variant.stamp: FORCE
	@if [ "$$(cat variant.stamp 2>/dev/null)" != "$(VARIANT)" ]; then echo $(VARIANT) > variant.stamp; fi

FORCE:

%.o: %.f90 $(LIB)/libufvm.a variant.stamp
	$(F90) $(F90FLAGS) -c $<

# the library depends on the sources of src: a changed source
# rebuilds and installs the library of this variant before the
# programs compile against it (COARRAY is read from the environment
# by src/Makefile through Makefile.in)
$(LIB)/libufvm.a: $(wildcard ../../src/*.f90) ../../src/OBJECTS
	$(MAKE) -C ../../src clean
	$(MAKE) -C ../../src
	$(MAKE) -C ../../src install

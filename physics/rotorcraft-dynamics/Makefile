include ../../Makefile.in

#=====================================================================#
# THE DYNAMICS OF A ROTORCRAFT: rigid bodies and joints stated as
# fields of time (multibody_dynamics), the triple pendulum of the
# dissertation, and the rotor. src/ is the mathematics and the
# computing; a program here states the bodies, the joints, the
# Lagrangian and calls minimize.
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

MODULES  = multibody_dynamics.o
PROGRAMS = triple_pendulum rotorcraft_dynamics

all: $(PROGRAMS)

$(PROGRAMS): %: %.o $(MODULES) $(LIB)/libufvm.a
	$(F90) $< $(MODULES) -o $@ $(LDFLAGS)

$(PROGRAMS:=.o): $(MODULES)

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

include ../../Makefile.in

#=====================================================================#
# THE TAYLOR-GREEN VORTEX on the periodic box, in two dimensions and in
# three. src/ is the mathematics and the computing; a program here
# states the fields, the laws, the conditions, the objective and the
# Lagrangian, and calls minimize.
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

PROGRAMS = taylor_green_vortex taylor_green_vortex_3d

all: $(PROGRAMS)

$(PROGRAMS): %: %.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

# the periodic boxes: box.geo meshed by gmsh; the box in three
# dimensions written by box_3d.py, since gmsh's transfinite volume
# does not accept a face copied from its periodic master
box.msh: box.geo
	gmsh -2 box.geo -o box.msh > gmsh.log

box_3d.msh: box_3d.py
	python3 box_3d.py 8

clean:
	@rm -f *.o *.mod *.so $(PROGRAMS) *.vtu gmsh.log variant.stamp

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

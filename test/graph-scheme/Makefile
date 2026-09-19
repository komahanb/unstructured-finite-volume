include ../../Makefile.in

#=====================================================================#
# The two programs on the continuous/discrete residual interface: the
# two-stage third-order SDIRK tableau and the Taylor-Green vortex;
# and the earlier SDIRK program on the expression and the
# residual_operator directly.
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

all: dirk_coefficients taylor_green_vortex taylor_green_study run

dirk_coefficients: dirk_coefficients.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

taylor_green_vortex: taylor_green_vortex.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

# the same statement under chosen parameters, convergence printed
taylor_green_study: taylor_green_study.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

run: test.o $(LIB)/libufvm.a
	$(F90) test.o -o $@ $(LDFLAGS)

# the periodic box, meshed by gmsh from its geometry
box.msh: box.geo
	gmsh -2 box.geo -o box.msh > gmsh.log

clean:
	@rm -f *.o *.mod *.so run dirk_coefficients taylor_green_vortex taylor_green_study *.vtu gmsh.log variant.stamp

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

include ../../Makefile.in

#=====================================================================#
# The studies, the audits and the derivations of schemes on the
# continuous/discrete residual interface; the problems they read are
# stated in physics/.
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

all: dirk_coefficients taylor_green_study van_der_pol_continuous_adjoint van_der_pol_study van_der_pol_design_audit euler_bernoulli_beam_grid run

dirk_coefficients: dirk_coefficients.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

van_der_pol_continuous_adjoint: van_der_pol_continuous_adjoint.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

van_der_pol_study: van_der_pol_study.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

# the expansion along the design audited: the derivatives of J by the
# state against the adjoint's, the Hessian read twice, the conditions
# bare; the program stops when a difference exceeds 1e-10
van_der_pol_design_audit: van_der_pol_design_audit.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

# the points of the beamline graded by Newton on the estimated error
euler_bernoulli_beam_grid: euler_bernoulli_beam_grid.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

# the same statement under chosen parameters, convergence printed
taylor_green_study: taylor_green_study.o $(LIB)/libufvm.a
	$(F90) $< -o $@ $(LDFLAGS)

run: test.o $(LIB)/libufvm.a
	$(F90) test.o -o $@ $(LDFLAGS)

# the periodic box of the study, meshed by gmsh from the geometry of
# the problem
box.msh: ../../physics/taylor_green/box.geo
	gmsh -2 ../../physics/taylor_green/box.geo -o box.msh > gmsh.log

clean:
	@rm -f *.o *.mod *.so run dirk_coefficients taylor_green_study van_der_pol_continuous_adjoint van_der_pol_study van_der_pol_design_audit euler_bernoulli_beam_grid *.vtu gmsh.log variant.stamp

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
# programs compile against it, so that a build here never reads a
# module interface older than the sources (COARRAY is read from the
# environment by src/Makefile through Makefile.in)
$(LIB)/libufvm.a: $(wildcard ../../src/*.f90) ../../src/OBJECTS
	$(MAKE) -C ../../src clean
	$(MAKE) -C ../../src
	$(MAKE) -C ../../src install

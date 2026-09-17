include ../../Makefile.in

#=====================================================================#
# The two programs on the continuous/discrete residual interface: the
# two-stage third-order SDIRK tableau and the Taylor-Green vortex;
# and the earlier SDIRK program on the expression and the
# residual_operator directly.
# -Wtarget-lifetime checks pointer target lifetimes at compile time.
#=====================================================================#

# Ufvm (-I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a

all: dirk_coefficients taylor_green_vortex run

dirk_coefficients: dirk_coefficients.o
	$(F90) $< -o $@ $(LDFLAGS)

taylor_green_vortex: taylor_green_vortex.o
	$(F90) $< -o $@ $(LDFLAGS)

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

# the periodic box, meshed by gmsh from its geometry
box.msh: box.geo
	gmsh -2 box.geo -o box.msh > gmsh.log

clean:
	@rm -f *.o *.mod *.so run dirk_coefficients taylor_green_vortex *.vtu gmsh.log

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

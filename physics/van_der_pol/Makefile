include ../../Makefile.in

#=====================================================================#
# Build the Van der Pol program with the mathematics and computation
# supplied by src/. The program states the fields, required relations
# and data, then the objective and Lagrangian. Calls to minimize use
# that construction to determine motion, adjoints, sensitivities along
# the design and the design itself.
# Use -Wtarget-lifetime to check pointer target lifetimes at compile
# time.
#=====================================================================#

# Use the modules and library for the chosen build: ../../lib for
# serial and lib_par for COARRAY=lib. Makefile.in supplies the coarray
# link flags; -I selects the corresponding .mod files.
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

# Recompile when the library changes because its module interfaces
# may have changed. Record the build variant as well, so switching
# between serial and coarray builds also recompiles the objects.
variant.stamp: FORCE
	@if [ "$$(cat variant.stamp 2>/dev/null)" != "$(VARIANT)" ]; then echo $(VARIANT) > variant.stamp; fi

FORCE:

%.o: %.f90 $(LIB)/libufvm.a variant.stamp
	$(F90) $(F90FLAGS) -c $<

# Rebuild and install this variant's library when its sources change,
# before compiling the program against it. src/Makefile reads COARRAY
# from the environment through Makefile.in.
$(LIB)/libufvm.a: $(wildcard ../../src/*.f90) ../../src/OBJECTS
	$(MAKE) -C ../../src clean
	$(MAKE) -C ../../src
	$(MAKE) -C ../../src install

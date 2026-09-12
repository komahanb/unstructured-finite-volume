include ../../Makefile.in

program = run

#=====================================================================#
# Object Files
#=====================================================================#

include OBJECTS

# Ufvm (-I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/
LDFLAGS += -L../../lib/ ../../lib/libufvm.a

#=====================================================================#
# Targets
#=====================================================================#

all: $(program) refusal restriction_refusal

$(program): $(objects)
	$(F90) $(objects) -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) refusal.o -o $@ $(LDFLAGS)

restriction_refusal: restriction_refusal.o
	$(F90) restriction_refusal.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so $(program) refusal restriction_refusal

#=====================================================================#
# Rule
#=====================================================================#

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

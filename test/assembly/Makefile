include ../../Makefile.in

program = run

#=====================================================================#
# Object Files
#=====================================================================#

include OBJECTS

# Ufvm
F90FLAGS += -I../../lib/
LDFLAGS += -L../../lib/ ../../lib/libufvm.a -llapack

#=====================================================================#
# Targets
#=====================================================================#

all: $(program) refusal

$(program): $(objects)
	$(F90) $(objects) -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) refusal.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so $(program) refusal

#=====================================================================#
# Rule
#=====================================================================#

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

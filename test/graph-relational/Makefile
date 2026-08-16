include ../../Makefile.in

#=====================================================================#
# The relational view over the production kernel.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: run lifetime

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

lifetime: lifetime.o
	$(F90) lifetime.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so run lifetime

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

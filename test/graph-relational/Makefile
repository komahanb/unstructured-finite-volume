include ../../Makefile.in

#=====================================================================#
# The relational view over the production kernel.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: run

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

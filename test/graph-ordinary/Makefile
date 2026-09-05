include ../../Makefile.in

#=====================================================================#
# One program: the law suite. -Wtarget-lifetime reports pointer
# targets going out of scope at compile time.
#=====================================================================#

# Ufvm (-I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -Wtarget-lifetime
LDFLAGS += -L../../lib/ ../../lib/libufvm.a

all: run

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)



clean:
	@rm -f *.o *.mod *.so run

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

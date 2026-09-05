include ../../Makefile.in

#=====================================================================#
# The partition / assembly law suite, for more than one part.
#=====================================================================#

# Ufvm (-I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/
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

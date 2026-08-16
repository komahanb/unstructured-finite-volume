include ../../Makefile.in

#=====================================================================#
# The identity maps: one storage law over set_map and inclusion_map.
# -Wtarget-lifetime keeps the borrows honest - there are none, and
# after the gate there are no TARGETs to be honest about either.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: lifetime

lifetime: lifetime.o
	$(F90) lifetime.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so lifetime

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

include ../../Makefile.in

#=====================================================================#
# Prototypes for the member-set and relation maps. Analysis only:
# nothing here is production, and nothing in src/ was changed to make
# these experiments convenient.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: set scale

set: set.o
	$(F90) set.o -o $@ $(LDFLAGS)

scale: scale.o
	$(F90) scale.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so set scale

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

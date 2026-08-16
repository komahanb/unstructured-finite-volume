include ../../Makefile.in

#=====================================================================#
# The declared embedding, over the kernel.
# -Wtarget-lifetime keeps the borrows honest - there are none.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: run refusal

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) refusal.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

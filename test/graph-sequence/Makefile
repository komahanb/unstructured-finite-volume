include ../../Makefile.in

#=====================================================================#
# The sequence view: the law suite, the refusals, and the two
# candidate API forms compared in candidates/.
#=====================================================================#

F90FLAGS += -I../../lib/ -Wtarget-lifetime
LDFLAGS  += -L../../lib/ ../../lib/libufvm.a

all: run refusal

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) refusal.o -o $@ $(LDFLAGS)

clean:
	@rm -f *.o *.mod *.so run refusal candidates/*.o candidates/*.mod candidates/compare

print-f90:
	@echo $(F90)

print-std:
	@echo $(FSTD)

.SUFFIXES: .f90 .o
.PHONY: all clean print-f90 print-std

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

include ../../Makefile.in

F90FLAGS += -I../../lib/
LDFLAGS += ../../lib/libufvm.a

all: run

run: test.o
	$(F90) test.o -o $@ $(LDFLAGS)

test.o: test.f90
	$(F90) $(F90FLAGS) -c $<

clean:
	rm -f *.o *.mod run

.PHONY: all clean

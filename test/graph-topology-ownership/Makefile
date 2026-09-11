include ../../Makefile.in

F90FLAGS += -I../../lib/ -I../../src/
LDFLAGS += ../../lib/libufvm.a

all: run traversal refusal

run: test.o
	$(F90) $< -o $@ $(LDFLAGS)

refusal: refusal.o
	$(F90) $< -o $@ $(LDFLAGS)

traversal: traversal.o allocation_counter.o
	$(F90) $^ -o $@ $(LDFLAGS) -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc

test.o: test.f90
	$(F90) $(F90FLAGS) -c $<

refusal.o: refusal.f90
	$(F90) $(F90FLAGS) -c $<

traversal.o: traversal.F90
	$(F90) $(F90FLAGS) -O3 -c $<

allocation_counter.o: allocation_counter.c
	$(CC) -O2 -c $<

refusals:
	python3 compile_refusals.py $(F90) $(F90FLAGS)

clean:
	rm -f *.o *.mod run traversal refusal

.PHONY: all clean refusals

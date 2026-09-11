include ../../Makefile.in
F90FLAGS += -I../../lib/
LDFLAGS += ../../lib/libufvm.a

all: run refusal

run: linear_rules.o test.o
	$(F90) $^ -o $@ $(LDFLAGS)

refusal: linear_rules.o refusal.o
	$(F90) $^ -o $@ $(LDFLAGS)

test.o refusal.o: linear_rules.o

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean:
	rm -f *.o *.mod run refusal

.PHONY: all clean

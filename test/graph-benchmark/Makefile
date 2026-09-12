include ../../Makefile.in

program = run

#=====================================================================#
# Object Files
#=====================================================================#

include OBJECTS

# Ufvm (-I../../lib for the .mod files)
F90FLAGS += -I../../lib/ -I../../src/ -O2
LDFLAGS += -L../../lib/ ../../lib/libufvm.a

# The scaling programs count allocation calls through wrapped
# allocators and read peak resident size from the kernel.
WRAP = -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc
scaling = schedule subset schur
# horizon links the packed application modules compiled into GTI_MODULE_DIR
GTI_MODULE_DIR ?= .

#=====================================================================#
# Targets
#=====================================================================#

all: $(program) $(scaling)

$(program): $(objects)
	$(F90) $(objects) -o $@ $(LDFLAGS)

schedule: schedule.o linear_rules.o measurement.o resource_counter.o
	$(F90) $^ -o $@ $(LDFLAGS) $(WRAP)

subset: subset.o measurement.o resource_counter.o
	$(F90) $^ -o $@ $(LDFLAGS) $(WRAP)

schur: schur.o measurement.o resource_counter.o
	$(F90) $^ -o $@ $(LDFLAGS) $(WRAP)

horizon: horizon.o measurement.o resource_counter.o
	$(F90) $^ $(GTI_MODULE_DIR)/modules.o -o $@ $(LDFLAGS) $(WRAP)

horizon.o: horizon.f90 measurement.o
	$(F90) $(F90FLAGS) -I$(GTI_MODULE_DIR) -c $<

schedule.o: schedule.f90 measurement.o linear_rules.o
subset.o: subset.f90 measurement.o
schur.o: schur.f90 measurement.o

linear_rules.o: ../graph-execution/linear_rules.f90
	$(F90) $(F90FLAGS) -c $< -o $@

resource_counter.o: resource_counter.c
	$(CC) -O2 -c $<

clean:
	@rm -f *.o *.mod *.so $(program) $(scaling) horizon

#=====================================================================#
# Rule
#=====================================================================#

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

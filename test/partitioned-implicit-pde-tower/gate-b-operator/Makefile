include ../../../Makefile.in

#=====================================================================#
# Gate B: the topology-consuming action. The test-local shifted
# adapter stores no graph and delegates all topology to the
# production differential operator; GMRES appears here for the
# first time in this tower, and its host is proved load-bearing
# rather than assumed so. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run refusal

run: partitioned_pde_assert.o shifted_laplacian_fixture.o test.o
	$(F90) partitioned_pde_assert.o shifted_laplacian_fixture.o test.o -o $@ $(LDFLAGS)

refusal: partitioned_pde_assert.o shifted_laplacian_fixture.o refusal.o
	$(F90) partitioned_pde_assert.o shifted_laplacian_fixture.o refusal.o -o $@ $(LDFLAGS)

test.o: partitioned_pde_assert.o shifted_laplacian_fixture.o
refusal.o: partitioned_pde_assert.o shifted_laplacian_fixture.o

shifted_laplacian_fixture.o: ../common/shifted_laplacian_fixture.f90
	$(F90) $(F90FLAGS) -c $< -o $@

partitioned_pde_assert.o: ../common/partitioned_pde_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

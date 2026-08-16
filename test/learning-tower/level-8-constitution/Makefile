include ../../../Makefile.in

#=====================================================================#
# Learning level 8: the meaning suite and the refusals that must
# die. The constitution is TEST-LOCAL - laws bind here, outside
# src/, and its evaluation semantics are the learning tower's own:
# laws execute INTO the computed slots, the residual is the value
# at the home L locates. Framework modules this level may see: the
# carriers, the relations, the algebra, the container, the profile,
# the algorithms and the field - no solver, no legacy host, no
# minimization. The import gate holds the door.
#=====================================================================#

F90FLAGS += -I../../../lib/ -I../../../src/
LDFLAGS += -L../../../lib/ ../../../lib/libufvm.a -llapack

all: run refusal

run: learning_assert.o learning_constitution_fixture.o test.o
	$(F90) learning_assert.o learning_constitution_fixture.o test.o -o $@ $(LDFLAGS)

refusal: learning_assert.o learning_constitution_fixture.o refusal.o
	$(F90) learning_assert.o learning_constitution_fixture.o refusal.o -o $@ $(LDFLAGS)

test.o: learning_assert.o learning_constitution_fixture.o
refusal.o: learning_assert.o learning_constitution_fixture.o
learning_constitution_fixture.o: learning_assert.o

learning_assert.o: ../common/learning_assert.f90
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	@rm -f *.o *.mod *.so run refusal

.SUFFIXES: .f90 .o
.PHONY: all clean

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

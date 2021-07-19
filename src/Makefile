include ../Makefile.in

#=====================================================================#
# Object Files
#=====================================================================#

include OBJECTS

#=====================================================================#
# Library
#=====================================================================#

default: $(objects) 
	ar rvs libufvm.a $(objects)
	@echo "----------- static library: libufvm.a created -----------"
	$(F90) $(F90FLAGS) -shared -o libufvm.so $(objects)
	@echo "----------- shared library: libufvm.so created -----------"

install:
	mv *.o ../lib
	mv *.a ../lib
	mv *.so ../lib
	mv *.mod ../lib

clean:
	@rm -f ../lib/*.o ../lib/*.mod ../lib/*.so ../lib/*.a *.mod *.o

#=====================================================================#
# Rule
#=====================================================================#

.SUFFIXES: .f90 .o
.PHONY: all clean 

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

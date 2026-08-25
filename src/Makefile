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
	mkdir -p $(LIBDIR)
	mv *.o $(LIBDIR)
	mv *.a $(LIBDIR)
	mv *.so $(LIBDIR)
	mv *.mod $(LIBDIR)

clean:
	@rm -f $(LIBDIR)/*.o $(LIBDIR)/*.mod $(LIBDIR)/*.so $(LIBDIR)/*.a *.mod *.o

#=====================================================================#
# Rule
#=====================================================================#

.SUFFIXES: .f90 .o
.PHONY: all clean 

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

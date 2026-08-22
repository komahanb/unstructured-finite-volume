# maths/: the technical report, its demonstrations and its experiments.
#
#   make              the report, maths/report.pdf
#   make demo         build lib/ then run demo/law (slice 1) and demo/law2 (slice 2)
#   make experiments  the time-order study and the Poisson study (needs gmsh for the meshes)
#   make oracle       restore the verification suites from ddcf592 and run them (maths/oracle.sh)
#   make clean

include ../Makefile.in

F90FLAGS += -I../lib/
LDFLAGS  += ../lib/libufvm.a

all: report.pdf

report.pdf: report.tex
	latexmk -pdf -interaction=nonstopmode report.tex

demo: ../lib/libufvm.a demo/law demo/law2
	cd demo && ./law && ./law2

experiments: ../lib/libufvm.a demo/time_orders demo/poisson demo/mesh/square-80.msh
	cd demo && ./time_orders && ./poisson

oracle:
	bash oracle.sh

../lib/libufvm.a:
	cd .. && ./build.sh

demo/law: demo/law.f90
	$(F90) $(F90FLAGS) demo/law.f90 -o demo/law $(LDFLAGS)

demo/law2: demo/law2.f90 demo/toy_forms.f90
	$(F90) $(F90FLAGS) -Jdemo -c demo/toy_forms.f90 -o demo/toy_forms.o
	$(F90) $(F90FLAGS) -Idemo demo/law2.f90 demo/toy_forms.o -o demo/law2 $(LDFLAGS)

demo/time_orders: demo/time_orders.f90
	$(F90) $(F90FLAGS) -O2 -Jdemo demo/time_orders.f90 -o demo/time_orders $(LDFLAGS)

demo/poisson: demo/poisson.f90
	$(F90) $(F90FLAGS) -O2 demo/poisson.f90 -o demo/poisson $(LDFLAGS)

demo/mesh/square-80.msh:
	mkdir -p demo/mesh
	for n in 10 20 40 80; do python3 ../meshgen/generate.py square-$$n demo/mesh/square-$$n.msh; done

clean:
	latexmk -C report.tex
	rm -f build.log demo/*.o demo/*.mod demo/law demo/law2 demo/time_orders demo/poisson
	rm -rf oracle

.PHONY: all demo experiments oracle clean

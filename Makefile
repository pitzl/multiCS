
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)

multics: multics.cc Makefile
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) -o multics multics.cc $(ROOTLIBS)
	@echo "done: multics (with ROOT)"

multicspp: multics.cpp Makefile
	g++ -O2 -Wall -Wextra -o multicspp multics.cpp
	@echo "done: multicspp (without ROOT)"

multicsf: multics.f90 Makefile
	gfortran -O2 -o multicsf multics.f90
	@echo "done: multicsf"

edgesc: edgesc.cc Makefile
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) -o edgesc edgesc.cc $(ROOTLIBS)
	@echo "done: edgesc (with ROOT)"

simloss: simloss.cc Makefile
	g++ -O2 -Wall -Wextra $(ROOTCFLAGS) -o simloss simloss.cc $(ROOTLIBS)
	@echo "done: simloss (with ROOT)"

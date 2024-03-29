ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS  += $(ROOTCFLAGS)
LIBS       = $(ROOTLIBS) 
GLIBS      = $(ROOTGLIBS)
GXX	   = /usr/bin/g++ -Wall -O3

all: RKnDemo

RKnTest:  RKnTest.cpp RKn.o
	$(GXX) -o RKnTest RKnTest.cpp  RKn.o $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)

RKnDemo:  RKnDemo.cpp RKn.o
	$(GXX) -o RKnDemo RKnDemo.cpp RKn.o $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS) 

RKn.o: RKn.cpp RKn.hpp
	$(GXX) -c -o RKn.o RKn.cpp $(ROOTCFLAGS) 

clean:
	rm -f RKnDemo RKnTest RKn.o *~

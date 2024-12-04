CXX=g++

LDLIBS=`root-config --libs`
CPPFLAGS=`root-config --cflags`

CPPFLAGS += -Wall -std=c++17 -O3

all: analysis
	@echo Compiling finished.

analysis: analysis.o cluster.o
	$(CXX) analysis.o cluster.o -o analysis $(CPPFLAGS) $(LDLIBS) 

analysis.o: analysis.cxx
	$(CXX) -c analysis.cxx -o analysis.o $(CPPFLAGS) $(LDLIBS) 

cluster.o: cluster.cxx cluster.hh
	$(CXX) -c cluster.cxx -o cluster.o

clean:
	$(RM) *.o *~ analysis

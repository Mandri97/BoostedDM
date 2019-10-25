CXX=g++

LDLIBS=`root-config --libs`
CPPFLAGS=`root-config --cflags`

CPPFLAGS += -Wall -std=c++11 -O3

all: analysis
	@echo Compiling finished.

analysis: analysis.o event.o cut.o
	$(CXX) analysis.o event.o cut.o -o analysis $(CPPFLAGS) $(LDLIBS) 

analysis.o: analysis.cxx
	$(CXX) -c analysis.cxx -o analysis.o $(CPPFLAGS) $(LDLIBS) 

cut.o: cut.cxx cut.hh
	$(CXX) -c cut.cxx -o cut.o

event.o: event.cxx event.hh
	$(CXX) -c event.cxx -o event.o


clean:
	$(RM) *.o *~ analysis

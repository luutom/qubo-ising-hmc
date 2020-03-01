## Comments
##

### GNU 
### change your compiler options here. . .
CXX=g++-9
CXXFLAGS=-std=c++11 -O2 -fopenmp -Wpedantic

### WARNINGS
CXXFLAGS+=-g -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter

### Lapack
### point to your local installation of lapack
CXXFLAGS+=-I/usr/local/Cellar/lapack/3.8.0_2/include -L/usr/local/Cellar/lapack/3.8.0_2/lib
LIBS=-llapacke -llapack

HDRS=quboIsingHMC.h
OBJS=quboIsingHMC.o

all: anneal

anneal: $(OBJS) anneal.o
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) anneal.o $(LIBS)

libquboIsingHMC.a:
	$(CXX) -w -c quboIsingHMC.cpp
	ar rcs libquboIsingHMC.a $(OBJS)

libisingHMC.so:
	$(CXX) -w -fPIC -c quboIsingHMC.cpp
	$(CXX) -shared -o libquboIsingHMC.so $(LIBS) $(OTHER_LIBS)

install: anneal

archive:
	tar -zcvf quboIsingHMC.tar.gz *.cpp *.h Makefile 

%.o: %.cc $(HDRS)
	$(CXX) -w $(CXXFLAGS) -c $<

clean:
	rm -rf *.o *~ 

distclean:
	rm -rf *.exe *.o *~ libquboIsingHMC.so libquboIsingHMC.a anneal



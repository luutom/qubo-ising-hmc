## Skeleton makefile
##

### GNU 
### change your compiler options here. . .
CXX=g++-10
CXXFLAGS=-std=c++11 -O2 -fopenmp -Wpedantic

### WARNINGS
CXXFLAGS+=-g -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter

HDRS=quboIsingHMC.h
OBJS=quboIsingHMC.o

all: anneal.exe

anneal.exe: $(OBJS) anneal.o
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) anneal.o $(LIBS)

libquboIsingHMC.a:
	$(CXX) -w -c quboIsingHMC.cpp
	ar rcs libquboIsingHMC.a $(OBJS)

libisingHMC.so:
	$(CXX) -w -fPIC -c quboIsingHMC.cpp
	$(CXX) -shared -o libquboIsingHMC.so $(LIBS) $(OTHER_LIBS)

install: anneal.exe

archive:
	tar -zcvf quboIsingHMC.tar.gz *.cpp *.h Makefile 

%.o: %.cc $(HDRS)
	$(CXX) -w $(CXXFLAGS) -c $<

clean:
	rm -rf *.o *~ 

distclean:
	rm -rf *.exe *.o *~ libquboIsingHMC.so libquboIsingHMC.a anneal



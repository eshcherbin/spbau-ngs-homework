CXX ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra -Wshadow -Werror --std=c++17

all: debrujin

debrujin: debrujin.o
	$(CXX) debrujin.o -o debrujin

debruijn.o: debrujin.cpp debrujin.hpp config.hpp
	$(CXX) debrujin.cpp $(CXXFLAGS) -c -o debrujin.o

graph.svg: graph.dot
	dot -Tsvg graph.dot -o graph.svg

clean:
	rm -f *.o debrujin edges.fasta graph.{dot,svg}

.PHONY: clean all graph.png

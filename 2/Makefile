CXX ?= g++
CXXFLAGS ?= -g -Wall -Wextra -Wshadow -Werror --std=c++17

all: debrujin

debrujin: debrujin.o
	$(CXX) debrujin.o -o debrujin

debruijn.o: debrujin.cpp debrujin.hpp utils.hpp
	$(CXX) debrujin.cpp $(CXXFLAGS) -c -o debrujin.o

graph.png: graph.dot
	dot -Tpng graph.dot -o graph.png

clean:
	rm -f *.o debrujin edges.fasta graph.{dot,png}

.PHONY: clean all graph.png

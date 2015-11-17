CXX=g++
CPPFLAGS=-g -Wall -Wextra $(shell root-config --cflags --glibs)

build: experiment analysis

experiment: detector.hpp experiment.cpp
	$(CXX) $(CPPFLAGS) experiment.cpp -o experiment

analysis: detector.hpp analysis.cpp
	$(CXX) $(CPPFLAGS) analysis.cpp -o analysis

clean:
	-rm -f experiment analysis

rebuild: clean build

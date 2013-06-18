
# author: Benjamin Aschenbrenner
# year: 2013
# Feel free to copy and modify just include this header.

CXX_COMPILER = g++
CXX_FLAGS = -std=c++11 -O2 -o

BUILD_PATH = bin/
EXECUTABLE_NAME = main

PLOTFILE = plotcommands.gp

#build main
all:
	mkdir -p bin
	tput reset
	$(CXX_COMPILER) $(CXX_FLAGS) $(BUILD_PATH)$(EXECUTABLE_NAME) magnetCalibration.cpp

test: all
	tput reset
	./bin/main
	gnuplot $(PLOTFILE)

run: run
	tput reset
	./bin/$(EXECUTABLE_NAME)
	gnuplot $(PLOTFILE)

clean:
	rm -rf bin


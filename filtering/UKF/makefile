
# author: Benjamin Aschenbrenner
# year: 2013
# Feel free to copy and modify just include this header.


CXX_COMPILER = g++
#CXX_FLAGS = -std=c++11 -O0 -g3 -ggdb -Wextra -o
CXX_FLAGS = -O3 -Wextra -o
LINK_FLAGS = -L /usr/lib/ -lboost_system -lboost_thread -lpthread

BUILD_PATH = bin/
EXECUTABLE_NAME = main

PLOTFILE = plotcommands.gp

#build main
all:
	mkdir -p bin
	tput reset
	$(CXX_COMPILER) test.cpp $(CXX_FLAGS) $(BUILD_PATH)$(EXECUTABLE_NAME) $(LINK_FLAGS) 

test: all
	tput reset
	./bin/main

run: 
	tput reset
	./bin/$(EXECUTABLE_NAME)
#	gnuplot $(PLOTFILE)

clean:
	rm -rf bin


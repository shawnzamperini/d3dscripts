CXX = g++
CFLAGS = -Wall -std=c++0x
DEPS = read_input.h
OBJ = kiss.o read_input.o

kiss: kiss.o read_input.o plasma_grid.o
	g++ -std=c++11 -Wall -o kiss kiss.o read_input.o plasma_grid.o

kiss.o: kiss.cpp
	g++ -std=c++11 -Wall -c kiss.cpp

read_input.o: read_input.cpp read_input.h
	g++ -std=c++11 -Wall -c read_input.cpp

plasma_grid.o: plasma_grid.cpp plasma_grid.h
	g++ -std=c++11 -Wall -c plasma_grid.cpp

all:
	kiss

clean:
	rm -f kiss kiss.o

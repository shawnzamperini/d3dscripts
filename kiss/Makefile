CXX = g++
CFLAGS = -Wall -std=c++0x
DEPS = read_input.h
OBJ = kiss.o read_input.o

kiss: kiss.o read_input.o plasma_grid.o follow_impurities.o
	g++ -g -std=c++11 -Wall -o kiss kiss.o read_input.o plasma_grid.o follow_impurities.o

kiss.o: kiss.cpp
	g++ -g -std=c++11 -Wall -c kiss.cpp

read_input.o: read_input.cpp read_input.h
	g++ -g -std=c++11 -Wall -c read_input.cpp

plasma_grid.o: plasma_grid.cpp plasma_grid.h
	g++ -g -std=c++11 -Wall -c plasma_grid.cpp

follow_impurities.o: follow_impurities.cpp follow_impurities.h
	g++ -g -std=c++11 -Wall -c follow_impurities.cpp

all:
	kiss

clean:
	rm -f kiss kiss.o

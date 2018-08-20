main: main.cpp
	g++ -O3 -I ../lib/include -L ../lib/lib -std=c++1y main.cpp kde/kde.c -o main -lhts -lgsl -lgslcblas -lm -pthread

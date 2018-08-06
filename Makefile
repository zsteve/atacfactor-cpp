main: main.cpp
	g++ -O3 -I ../lib/include -L ../lib/lib -std=c++1y main.cpp -o main -lhts

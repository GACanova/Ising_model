output: main.o
	g++ -std=c++11 -O3 main.o -o out 

main.o: main.cpp
	g++ -c main.cpp 

clean:
	rm *.o out

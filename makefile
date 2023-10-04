
test.exe: DecisionT.o algorithm.o test.o
	g++ -g -std=c++11 -o test.exe DecisionT.o algorithm.o test.o


DecisionT.o : DecisionT.cpp 
	g++ -g -std=c++11 -c DecisionT.cpp

algorithm.o : algorithm.cpp
	g++ -g -std=c++11 -c algorithm.cpp

test.o: test.cpp
	g++ -g -std=c++11 -c test.cpp

clean:
	del *.o test.exe


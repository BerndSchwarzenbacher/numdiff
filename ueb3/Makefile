CXX = g++
FLAGS = -std=c++11 ../ode/bla/bla.so -O2
OBJECTS = ex12 ex13 ex14

all: $(OBJECTS)

%: %.cpp
	g++ $? $(FLAGS) -o $@

clean:
	rm $(OBJECTS)


CXX = g++
FLAGS = -std=c++11 ../ode/bla/bla.so
OBJECTS = ex12.o ex14.o

all: $(OBJECTS)

%.o: %.cpp
	g++ $? $(FLAGS) -o $@

clean:
	rm $(OBJECTS)


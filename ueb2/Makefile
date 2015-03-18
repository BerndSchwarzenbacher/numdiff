CXX = g++
FLAGS = -std=c++11 ../ode/bla/bla.so
OBJECTS = ex06.o ex08.o

all: $(OBJECTS)

%.o: %.cpp
	g++ $? $(FLAGS) -o $@

clean:
	rm $(OBJECTS)


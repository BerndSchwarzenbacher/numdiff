CXX = g++
FLAGS = -std=c++11 -c -fPIC
OBJECTS = exception.o calcinverse.o localheap.o

all: bla.so

bla.so: $(OBJECTS)
	$(CXX) -shared $? -o $@

%.o: %.cpp
	g++ $? $(FLAGS) -o $@

clean:
	rm $(OBJECTS) bla.so


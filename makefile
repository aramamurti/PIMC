CXX = mpic++
CXX_FLAGS = -std=c++14 -Ofast
LD_FLAGS = -lgsl -lgslcblas

EXEC = PIMC
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LD_FLAGS) -o $(EXEC)

%.o:%.cpp
	$(CXX) -c $(CXX_FLAGS) $< -o $@

clean:
	rm -f $(OBJECTS)


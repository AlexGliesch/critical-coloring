TARGET = critcol
GENERATOR = generateimperfectgraph

override CFLAGS +=-Wall -Wextra -Wfatal-errors -std=c++17 -MD -MP -O3 -fconcepts
CPP = g++
LIBS = -lpthread -lstdc++ -lboost_program_options -lboost_system -lboost_filesystem 
SRC = $(filter-out $(GENERATOR).cpp, $(wildcard *.cpp btdsatur/*.cpp hybridea/*.cpp cliques/mntshao.cpp))
HEADERS = $(wildcard *.h *.inl)
OBJS = $(SRC:.cpp=.o)

all: $(TARGET) $(GENERATOR)

release: CFLAGS += -DNDEBUG
release: $(TARGET) $(GENERATOR)

pre.h.gch: pre.h random.h
	$(CPP) -O3 -fconcepts -std=c++17 pre.h 
	
$(TARGET): $(OBJS)
	$(CPP) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIBS)

$(GENERATOR): $(GENERATOR).o
	$(CPP) $(CFLAGS) $(INCLUDES) -o $(GENERATOR) $(GENERATOR).o btdsatur/*.o hybridea/*.o cliques/*.o $(LIBS)		

-include $(SRC:.cpp=.d)

%.o: %.cpp pre.h.gch
	$(CPP) $(CFLAGS) -include pre.h -c $< -o $@

clean:
	rm -f *.o *.d btdsatur/*.o btdsatur/*.d hybridea/*.o hybridea/*.d cliques/*.d cliques/*.o $(TARGET)  $(GENERATOR) *.stackdump

cleanall:
	$(MAKE) clean && rm pre.h.gch
TARGET = ibr

override CFLAGS +=-Wall -Wextra -Wfatal-errors -std=gnu++1z -MD -MP -O3
CPP = g++
LIBS = -lpthread -lstdc++ -lboost_program_options -lboost_system -lboost_filesystem
SRC = $(wildcard *.cpp ../../btdsatur/*.cpp ../../hybridea/*.cpp ../../cliques/mntshao.cpp)
HEADERS = $(wildcard *.h *.inl)
OBJS = $(SRC:.cpp=.o)

all: $(TARGET) 

release: CFLAGS += -DNDEBUG
release: $(TARGET) 

pre.h.gch: ../../pre.h ../../random.h
	$(CPP) -O3 -std=gnu++1z ../../pre.h -o pre.h.gch
	
$(TARGET): $(OBJS)
	$(CPP) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) $(LIBS)

-include $(SRC:.cpp=.d)

%.o: %.cpp pre.h.gch
	$(CPP) $(CFLAGS) -include ../../pre.h -c $< -o $@

clean:
	rm -f *.o *.d ../../btdsatur/*.o ../../btdsatur/*.d ../../hybridea/*.o ../../hybridea/*.d $(TARGET) *.stackdump

cleanall:
	$(MAKE) clean && rm pre.h.gch
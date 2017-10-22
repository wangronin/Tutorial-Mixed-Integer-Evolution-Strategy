FLAGS = -std=c++0x -O3 -DNDEBUG -Wall
CC=g++
OBJS=$(notdir $(SOURCES:.cpp=.o))
SOURCES=$(wildcard *.cpp) 

optical: ${OBJS}
	$(CC) $^ -o $@ 

%.o: %.cpp
	$(CC) $(FLAGS) -c $< -o $@

clean:
	rm -f ${OBJS} optical

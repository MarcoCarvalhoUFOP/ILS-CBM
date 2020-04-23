# Macros:
CPP= g++
CFLAGS= -O3 -std=c++11
MPFLAGS= -O3 -std=c++11
LFLAGS= -lm
MAIN= testManager.cpp

# Targets:
all: compile

compile:
	$(CPP) $(CFLAGS) $(MAIN) -o run

run:
	./run -pp 1 -s 3 -i 3 -alpha 0.2 -beta 0.5 -gamma 1 -p 1 -ls 3 -mh 1  -it 100 -r 1 -v 1 

# Remove:
clean:
	rm *.o *.gch run

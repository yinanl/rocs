CC = g++ -g -Wall
CFLAGS = -I.. -D_DEBUG

SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp, %.o, $(SRC))

.cpp.o:
	$(CC) $(CFLAGS) -c $<

all: $(testAll)

testAll: $(OBJ)
	$(CC) -o $@ $(OBJ)

exptad1: ExampleTAD1.o
	$(CC) -o $@ ExampleTAD1.o

exptad2: ExampleTAD2.o
	$(CC) -o $@ ExampleTAD2.o

exptad3: ExampleTAD3.o
	$(CC) -o $@ ExampleTAD3.o

disp:
	@echo $(OBJ)

clean:
	-rm *.o testAll

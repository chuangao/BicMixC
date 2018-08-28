CC=g++
unamestr := $(shell uname)
ifeq ($(unamestr),Linux)
	CFLAGS=-g -O2 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math -fopenmp
endif
ifeq ($(unamestr),Darwin)
	CFLAGS=-g -O2 -DNDEBUG -fpic -pipe -march=x86-64 -ffast-math
endif
#CFLAGS=-O2 -march=x86-64 -ffast-math -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -g -pg -fopenmp
LINKER=-L${CURDIR}/lib/ -Wl,-rpath,${CURDIR}/lib/

INCLUDE=-isystem${CURDIR}/ -isystem${CURDIR}/include/

#LIBS=-lstdc++ -lgfortran -lgsl -lgslcblas -lpthread -lm -static
LIBS=-lgsl -lgslcblas -lm 

OBJECTS= BicMix.o

all: install myHeader.o BicMix clean

BicMix: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o BicMix $(LINKER) $(LIBS)
BicMix.o:  BicMix.cpp myHeader.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  BicMix.cpp
myHeader.o: myHeader.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c  myHeader.cpp
install:
	./install.sh
clean:
	rm *.o

CC=g++

SRCS:=$(shell ls *.cpp *.c)
OBJS:=$(SRCS:.cpp=.o)
OBJS:=$(OBJS:.c=.o)

CFLAGS=-O3 -std=c++0x -Wall 
#LFLAGS= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

OUT=a

all:$(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -openmp -o $(OUT) 

.cpp.o:
	$(CC) $(CFLAGS) -c $< -fpermissive 

clean:
	rm $(OBJS) $(OUT)

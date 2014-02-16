# Sequential makefile for BC
CC = icc 
CFLAGS = 

AR = ar
ARFLAGS = cr
RANLIB = ranlib

LIB = -lm
INC = 
TARGET = bc

OBJS = bc.o utils.o generateGraphs.o betweennessCentrality.o 

.cpp.o: defs.h Makefile
	$(CC) $(INC) $(CFLAGS) -c $<


all: $(OBJS) defs.h Makefile
	$(CC) $(INC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIB)

clean: 
	rm -f *.o *~ $(TARGET) core*

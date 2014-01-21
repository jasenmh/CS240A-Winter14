CC = mpicc
CFLAGS =
LIBS = -lm
OBJECTS = hello_mpi.o
EXECUTABLE = hello_mpi

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o $(EXECUTABLE)

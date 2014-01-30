# Makefile for Cilk Plus inner products.  See source code for details.

CILKPP	= icc
LIBARG	= -O2
TARGET	= innerproduct 
SRC	= $(addsuffix .cpp,$(TARGET))

all: $(TARGET)

$(TARGET): $(SRC)
	$(CILKPP) $(SRC) $(LIBARG) -o $@

clean:
	rm -f $(TARGET)

# Makefile for standalone test of hrstrans
# R. Michaels, Oct 2017


# Choose the compiler.
GCC=g++
GLD=g++


ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
INCLUDES      = -I$(ROOTSYS)/include
CXX           = $(GCC)
CXXFLAGS      = -fno-exceptions -fpermissive -std=c++11 -fPIC $(INCLUDES)
LD            = $(GLD)
SOFLAGS       = -std=c++11  -shared 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 
LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS) /usr/lib64/libg2c.so.0

MAKEDEPEND    = $(GCC)

ALL_LIBS = $(LIBS) 

ifdef OPTIMIZE
   CXXFLAGS += -O
else
   CXXFLAGS += -g -ggdb
endif

SRC = THRSTrans.C

DEPS = $(SRC:.C=.d)
DEP  = $(SRC:.C=.d)
HEAD = $(SRC:.C=.h) 
OBJS = $(SRC:.C=.o)

install: all

all: libhrstrans.so 

libhrstrans.so: $(OBJS) $(HEAD)
	$(CXX) $(SOFLAGS) -O -o libhrstrans.so $(OBJS) $(ALL_LIBS)

clean:
	rm -f *.o core libhrstrans.so 

realclean:  clean
	rm -f *.tar  *~


%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<


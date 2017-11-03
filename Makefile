# Makefile for standalone test of hrstrans
# R. Michaels, Oct 2017


# Choose the compiler.
GCC=g++
GLD=g++


ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
INCLUDES      = -I$(ROOTSYS)/include
CXX           = $(GCC)
CXXFLAGS      = -Wall -fno-exceptions -fPIC $(INCLUDES)
LD            = $(GLD)
LDFLAGS       = 
SOFLAGS       = -shared 
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

all: xhrst libhrstrans.so 

xhrst: $(OBJS) $(SRC) $(HEAD) main.C 
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ main.C $(OBJS) $(ALL_LIBS)

libhrstrans.so: $(OBJS) $(HEAD)
	$(CXX) $(SOFLAGS) -O -o libhrstrans.so $(OBJS) $(ALL_LIBS)

clean:
	rm -f *.o core *Dict* $(PROGS) libhrstrans.so 

realclean:  clean
	rm -f *.d *.tar  *~


%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

%.d:	%.C 
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(MAKEDEPEND) -MM $(INCLUDES) -c $< \
                | sed '\''s%^.*\.o%$*\.o%g'\'' \
                | sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
                [ -s $@ ] || rm -f $@'

-include $(DEPS)

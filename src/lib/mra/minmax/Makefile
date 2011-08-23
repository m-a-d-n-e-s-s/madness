INCLUDES = -I/home/rjh/install/include -I.
LIBS = -L/home/rjh/install/lib -lqd -lm
CXXOPT = -g -Wall -Wno-deprecated -ansi -O0

#-ffast-math ... this breaks qd_real

CXXFLAGS = $(INCLUDES) $(CXXOPT)
CXX = g++

ex1:	ex1.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

yacobi:	yacobi.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

y:	y.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

yx:	yx.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)


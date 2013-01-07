INCLUDES = -I/home/rjh/install/include -I.
LIBS = -L/home/rjh/install/lib -lqd -lm
CXXOPT = -g -O3 -march=native
# -Wall -Wno-strict-aliasing -Wno-deprecated  -ffast-math -std=c++0x 
#CXXOPT = -mavx -g -O3 -Wall -Wno-strict-aliasing -Wno-deprecated  -ffast-math -std=c++0x -march=native
#CXXOPT = -g -Wall -Wno-deprecated -ansi -O0

#-ffast-math ... this breaks qd_real

CXXFLAGS = $(INCLUDES) $(CXXOPT)
CXX = g++

y_signed:	y_signed.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

y_varx:	y_varx.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

ex1:	ex1.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

yacobi:	yacobi.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

y:	y.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

yx:	yx.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)


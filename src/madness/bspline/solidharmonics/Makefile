INCLUDES = -I.
CXXFLAGS = -g -O3 -Wall $(INCLUDES)
LDLIBS = -lqd

solidharmonics:	solidharmonics.cc myqd.h gauleg.h
	$(CXX) $(CXXFLAGS) solidharmonics.cc -o $@ $(LDLIBS)

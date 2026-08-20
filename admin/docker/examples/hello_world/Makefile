# Makefile for building a MADNESS application against the installed container libraries
CXX ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall
CPPFLAGS ?= -I/usr/local/include
LDFLAGS ?= -L/usr/local/lib
LIBS ?= -lmadness -lxc -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

TARGET = hello_madness
SRCS = hello.cc
OBJS = $(SRCS:.cc=.o)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $^ $(LDFLAGS) $(LIBS) -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJS)

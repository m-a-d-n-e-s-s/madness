# Makefile for building a MADNESS application against the installed container
# libraries. CMakeLists.txt in this directory is the supported route -- it picks
# up MADNESS's public compile definitions and BLAS include paths from
# madness-config.cmake. This file spells the amd64/MKL link line out by hand;
# on a non-x86-64 image replace the mkl_* libraries with `-lopenblas -llapacke`.
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

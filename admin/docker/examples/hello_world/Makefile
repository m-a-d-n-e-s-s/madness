# Makefile for building a MADNESS application against the installed container
# libraries. CMakeLists.txt in this directory is the supported route -- it picks
# up MADNESS's public compile definitions and BLAS include paths from
# madness-config.cmake. This file spells the link line out by hand:
# on amd64 it links Intel MKL; on arm64 it links ARMPL.
ifeq ($(origin CXX),default)
  CXX = mpicxx
endif
CXX ?= mpicxx
CXXFLAGS ?= -std=c++20 -O2 -Wall
UNAME_M := $(shell uname -m)
ifneq (,$(filter $(UNAME_M),aarch64 arm64))
  CPPFLAGS ?= -I/usr/local/include -I/opt/arm/arm-performance-libraries/include
  LDFLAGS ?= -L/usr/local/lib -L/opt/arm/arm-performance-libraries/lib
  LIBS ?= -lmadness -lxsmm -lxc -ltbb -ltcmalloc_minimal -larmpl_lp64 -lamath -lastring -lpthread -lm -ldl
else
  CPPFLAGS ?= -I/usr/local/include -I/usr/include/mkl
  LDFLAGS ?= -L/usr/local/lib
  LIBS ?= -lmadness -lxsmm -lxc -ltbb -ltcmalloc_minimal -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
endif

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

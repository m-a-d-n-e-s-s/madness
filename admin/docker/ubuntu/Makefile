IMAGE ?= madness:latest
DOCKERFILE ?= Dockerfile
CONTEXT ?= ../../..

.PHONY: all build dev clean

all: build

build:
	docker build -t $(IMAGE) -f $(DOCKERFILE) $(CONTEXT)

dev:
	docker run -it --rm -v $$(pwd):/work -w /work $(IMAGE) /bin/bash

test:
	docker run --rm $(IMAGE) madqc --help
	docker run --rm $(IMAGE) bash -c "mpirun -np 2 -x MAD_NUM_THREADS=2 madqc --geometry=he"
	docker run --rm -v $$(pwd)/../examples/hello_world:/work -w /work --user "$$(id -u):$$(id -g)" $(IMAGE) bash -c "make clean && make && ./hello_madness && make clean"
	docker run --rm -v $$(pwd)/../examples/hello_world:/work -w /work --user "$$(id -u):$$(id -g)" $(IMAGE) bash -c "rm -rf build && cmake -B build -S . && cmake --build build && ./build/hello_madness && rm -rf build"

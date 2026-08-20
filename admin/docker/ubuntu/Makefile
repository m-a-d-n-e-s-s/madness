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
	docker run --rm $(IMAGE) madqc --version || true

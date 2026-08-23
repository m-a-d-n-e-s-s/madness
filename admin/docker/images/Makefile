repo ?= rjharrison
latest ?= 22.04

HOST_ARCH := $(shell uname -m)
ifeq ($(HOST_ARCH),aarch64)
  ARCH ?= arm64
else ifeq ($(HOST_ARCH),arm64)
  ARCH ?= arm64
else ifeq ($(HOST_ARCH),x86_64)
  ARCH ?= x86_64
else
  ARCH ?= $(HOST_ARCH)
endif

# Build configuration options
BUILD_TYPE ?= Release
BUILD_PARALLEL_JOBS ?= 2
CMAKE_EXTRA_ARGS ?=
MADNESS_REPO ?= https://github.com/m-a-d-n-e-s-s/madness
DOCKER_BUILD_FLAGS ?= --network=host

BUILD_ARGS = $(DOCKER_BUILD_FLAGS) \
             --build-arg ubuntuImage=ubuntu:$* \
             --build-arg BUILD_TYPE=$(BUILD_TYPE) \
             --build-arg BUILD_PARALLEL_JOBS=$(BUILD_PARALLEL_JOBS) \
             --build-arg CMAKE_EXTRA_ARGS="$(CMAKE_EXTRA_ARGS)" \
             --build-arg MADNESS_REPO="$(MADNESS_REPO)"

# Architecture-specific image tags
ALL_ARM64 = ${repo}/ubuntu-arm64\:22.04
ALL_X86_64 = ${repo}/ubuntu-x86_64\:22.04
ALL = ${repo}/ubuntu\:${latest}

# Build rules for arm64
${repo}/ubuntu-arm64\:%:
	docker build $(BUILD_ARGS) --build-arg ARCH=arm64 -f ../ubuntu/Dockerfile.arm64 ../ubuntu/ -t $@

# Build rules for x86_64
${repo}/ubuntu-x86_64\:%:
	docker build $(BUILD_ARGS) --build-arg ARCH=x86_64 -f ../ubuntu/Dockerfile.x86_64 ../ubuntu/ -t $@

# Unified / default build rule
${repo}/ubuntu\:%:
	docker build $(BUILD_ARGS) --build-arg ARCH=$(ARCH) -f ../ubuntu/Dockerfile ../ubuntu/ -t $@

# Convenience aliases for latest
${repo}/ubuntu-arm64: ${repo}/ubuntu-arm64\:${latest}
	docker tag ${repo}/ubuntu-arm64:${latest} $@

${repo}/ubuntu-x86_64: ${repo}/ubuntu-x86_64\:${latest}
	docker tag ${repo}/ubuntu-x86_64:${latest} $@

${repo}/ubuntu: ${repo}/ubuntu\:${latest}
	docker tag ${repo}/ubuntu:${latest} $@

arm64: ${repo}/ubuntu-arm64\:${latest}
x86_64: ${repo}/ubuntu-x86_64\:${latest}

all: arm64 x86_64

# Tar export targets
${repo}/ubuntu-arm64\:%.tar:
	mkdir -p ${repo}
	DOCKER_BUILDKIT=1 docker build $(BUILD_ARGS) --build-arg ARCH=arm64 -f ../ubuntu/Dockerfile.arm64 ../ubuntu/ -o - > $@

${repo}/ubuntu-x86_64\:%.tar:
	mkdir -p ${repo}
	DOCKER_BUILDKIT=1 docker build $(BUILD_ARGS) --build-arg ARCH=x86_64 -f ../ubuntu/Dockerfile.x86_64 ../ubuntu/ -o - > $@

${repo}/ubuntu\:%.tar:
	mkdir -p ${repo}
	DOCKER_BUILDKIT=1 docker build $(BUILD_ARGS) --build-arg ARCH=$(ARCH) -f ../ubuntu/Dockerfile ../ubuntu/ -o - > $@

all/tar: $(ALL_ARM64:%=%.tar) $(ALL_X86_64:%=%.tar)

# Push targets
push/arm64/%: ${repo}/ubuntu-arm64\:%
	docker push $?

push/x86_64/%: ${repo}/ubuntu-x86_64\:%
	docker push $?

push/%: ${repo}/ubuntu\:%
	docker push $?

push/all: push/arm64/22.04 push/x86_64/22.04

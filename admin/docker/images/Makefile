# Build and publish the MADNESS container images.
#
# The Dockerfile builds MADNESS from the surrounding source tree, so the build
# context is the repository root (the directory holding the top-level
# CMakeLists.txt), not admin/docker/ubuntu.
#
#   make all         # build an image per $(versions), tag the latest one
#   make build/24.04 # build one image
#   make all/tar     # export each image's filesystem as a tarball
#   make push/all    # build, tag and push everything
#
# N.B. targets are spelled `build/<version>` rather than `<repo>/ubuntu:<version>`
# because escaping the `:` in a make target name is not portable across make
# versions -- the previous `${repo}/ubuntu\:%` spelling never matched.

repo ?= rjharrison
latest ?= 24.04
versions ?= 24.04

dockerfile = ../ubuntu/Dockerfile
context = ../../..

# N.B. build/%, tar/% and push/% are deliberately *not* .PHONY: GNU make skips
# implicit- and pattern-rule search for phony targets, so listing them there
# would make `make all` silently do nothing. They never create a file of their
# own name, so they re-run every time anyway.
.PHONY: all all/tar push/all tag/latest

all: $(versions:%=build/%) tag/latest

build/%:
	docker build --build-arg UBUNTU_VERSION=$* -f $(dockerfile) $(context) -t $(repo)/ubuntu:$*

tar/%:
	mkdir -p $(repo)
	DOCKER_BUILDKIT=1 docker build --build-arg UBUNTU_VERSION=$* -f $(dockerfile) $(context) -o - > $(repo)/ubuntu-$*.tar

all/tar: $(versions:%=tar/%)

tag/latest: build/$(latest)
	docker tag $(repo)/ubuntu:$(latest) $(repo)/ubuntu

push/%: build/%
	docker push $(repo)/ubuntu:$*

push/all: $(versions:%=push/%) tag/latest
	docker push $(repo)/ubuntu

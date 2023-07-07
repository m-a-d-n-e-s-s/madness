repo = rjharrison
latest ?= 22.04
ALL = ${repo}/ubuntu\:22.04

${repo}/ubuntu\:%: build = docker build --build-arg ubuntuImage=ubuntu:$* -f ../ubuntu/Dockerfile ../ubuntu/

${repo}/ubuntu: ${repo}/ubuntu\:${latest}
	docker tag ${repo}/ubuntu:${latest} $@

${repo}/ubuntu\:%:
	${build} -t $@

${repo}/ubuntu\:%.tar:
	mkdir -p ${repo}
	DOCKER_BUILDKIT=1 ${build} -o - > $@

all: ${ALL} ${repo}/ubuntu

all/tar: $(ALL:%=%.tar)

push/%: ${repo}/ubuntu\:%
	docker push $?

push/all: push/22.04

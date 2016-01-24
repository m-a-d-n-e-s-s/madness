CMAKE := cmake-3.4.2
URL := http://www.cmake.org/files/v3.4/$(CMAKE).tar.gz
PREFIX := $(PWD)

all: cmake

wget := $(shell wget 2>/dev/null)
curl := $(shell curl 2>/dev/null)

ifneq  ($(wget), )
download := wget $(URL) -O $(CMAKE).tar.gz
else
download := curl $(URL) -o $(CMAKE).tar.gz
endif


$(CMAKE).tar.gz:
	$(download)

.PHONY: cmake
cmake: $(CMAKE).tar.gz
	tar -xzf $(CMAKE).tar.gz
	cd ./$(CMAKE)
	cd $(CMAKE) && ./configure --prefix=$(PREFIX) # --system-curl
	cd $(CMAKE) && make install



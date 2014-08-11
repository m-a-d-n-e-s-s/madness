SUBDIRS = gtest lib apps examples

libraries:
	for dir in lib apps; do $(MAKE) -C $$dir $@ ; done

install-libraries:
	for dir in lib ; do $(MAKE) -C $$dir $@ ; done

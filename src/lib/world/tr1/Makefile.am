include $(top_srcdir)/config/MakeGlobal.am


thisincludedir = $(includedir)/world/tr1
thisinclude_HEADERS = memory.h functional.h type_traits.h array.h 

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

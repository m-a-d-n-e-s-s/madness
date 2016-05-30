include $(top_srcdir)/config/MakeGlobal.am

if CAN_BUILD_DOC
.PHONY: docs
docs:
	$(DOXYGEN) doxygen.cfg


.PHONY: clean
clean::
	/bin/rm -rf html latex
	/bin/rm -rf *~

veryclean: clean

distclean: clean
	/bin/rm doxygen.cfg

install::
	-mkdir -p $(installroot)$(prefix)/doc
	-cp -r html $(installroot)$(prefix)/doc/

install_devel::

else
docs:
	doxygen is not available ... cannot build documentation ... sorry
endif


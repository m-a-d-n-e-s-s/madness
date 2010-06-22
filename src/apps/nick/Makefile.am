include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS     = projPsi test hello lineplot
projPsi_SOURCES  = projPsi.cc wavef.cc wavef.h hyp.cc hyp.h mpreal.cc
lineplot_SOURCES = lineplot.cc

test_SOURCES     = test.cc hyp.cc mpreal.cc
hello_SOURCES    = hello.cc wavef.cc wavef.h
projPsi_LDADD    = $(MRALIBS)
test_LDADD       = $(MRALIBS)
hello_LDADD      = $(MRALIBS)
lineplot_LDADD   = $(MRALIBS)

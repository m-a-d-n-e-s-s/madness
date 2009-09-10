include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS    = projPsi test hello
projPsi_SOURCES = projPsi.cc wavef.cc wavef.h hyp.cc hyp.h mpreal.cc

test_SOURCES    = test.cc hyp.cc mpreal.cc
hello_SOURCES   = hello.cc wavef.cc wavef.h
projPsi_LDADD   = $(MRALIBS)
test_LDADD      = $(MRALIBS)
hello_LDADD     = $(MRALIBS)

include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS    = projPsi test hello
projPsi_SOURCES = projPsi.cc wavef.cc coulcc36-f90.f wavef.h hyp.cc mpreal.cc
test_SOURCES    = test.cc coulcc36-f90.f hyp.cc mpreal.cc
hello_SOURCES   = hello.cc wavef.cc coulcc36-f90.f wavef.h
projPsi_LDADD   = $(MRALIBS) -lmpfr -lgmp
test_LDADD      = $(MRALIBS) -lmpfr -lgmp
hello_LDADD     = $(MRALIBS)

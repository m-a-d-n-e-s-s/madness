include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = projPsi test hello
projPsi_SOURCES = projPsi.cc wavef.cc coulcc36-f90.f wavef.h
test_SOURCES =  test.cc
hello_SOURCES =  hello.cc wavef.cc coulcc36-f90.f wavef.h
projPsi_LDADD = $(MRALIBS) 
test_LDADD = $(MRALIBS)
hello_LDADD = $(MRALIBS)
F77=gfortran

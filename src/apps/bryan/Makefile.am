include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

bin_PROGRAMS = driver driver2 test_ops

driver_SOURCES = driver.cc TDA.cc 
driver_LDADD = $(LIBCHEM) $(MRALIBS)


driver2_SOURCES = driver2.cc TDA2.cc 
driver2_LDADD = $(LIBCHEM) $(MRALIBS)

test_ops_SOURCES = test_ops.cc
test_ops_LDADD = $(LIBCHEM) $(MRALIBS)

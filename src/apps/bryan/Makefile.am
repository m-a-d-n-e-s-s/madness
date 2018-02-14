include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt


bin_PROGRAMS = driver driver2 driver3 driver4 test_ops


driver_SOURCES = driver.cc TDA.cc 
driver_LDADD = $(LIBCHEM) $(MRALIBS)

driver2_SOURCES = driver2.cc TDA2.cc
driver2_LDADD = $(LIBCHEM) $(MRALIBS)

driver3_SOURCES = driver3.cc TDA3.cc
driver3_LDADD = $(LIBCHEM) $(MRALIBS)

driver4_SOURCES = driver4.cc TDHF.cc
driver4_LDADD = $(LIBCHEM) $(MRALIBS)

test_ops_SOURCES = test_ops.cc
test_ops_LDADD = $(LIBCHEM) $(MRALIBS)



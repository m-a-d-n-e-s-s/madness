include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt


bin_PROGRAMS = driver driver2 driver3 driver4 driver5 test_ops test_ResponseFunction test_Gamma


driver_SOURCES = driver.cc TDA.cc 
driver_LDADD = $(LIBCHEM) $(MRALIBS)

driver2_SOURCES = driver2.cc TDA2.cc
driver2_LDADD = $(LIBCHEM) $(MRALIBS)

driver3_SOURCES = driver3.cc TDA3.cc
driver3_LDADD = $(LIBCHEM) $(MRALIBS)

driver4_SOURCES = driver4.cc TDHF.cc ../chem/polynomial.cc ../chem/gaussian.cc ../chem/NWChem.cc 
driver4_LDADD = $(LIBCHEM) $(MRALIBS)

driver5_SOURCES = driver5.cc TDHF2.cc ../chem/polynomial.cc ../chem/gaussian.cc ../chem/NWChem.cc
driver5_LDADD = $(LIBCHEM) $(MRALIBS)

test_ops_SOURCES = test_ops.cc
test_ops_LDADD = $(LIBCHEM) $(MRALIBS)

test_ResponseFunction_SOURCES = test_ResponseFunction.cc
test_ResponseFunction_LDADD = $(LIBCHEM) $(MRALIBS)

test_Gamma_SOURCES = test_Gamma.cc
test_Gamma_LDADD = $(LIBCHEM) $(MRALIBS)

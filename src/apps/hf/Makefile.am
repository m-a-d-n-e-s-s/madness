include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test_hf test_lattice test_coulomb test_hydro
lib_LIBRARIES = libHartreeFock.a

libHartreeFock_a_SOURCES = eigsolver.cc dft.cc dft.h  eigsolver.h poperator.h util.h

test_hf_SOURCES = test_hf.cc

test_hydro_SOURCES = test_hydro.cc

test_lattice_SOURCES = test_lattice.cc

test_coulomb_SOURCES = test_coulomb.cc

LDADD = libHartreeFock.a $(MRALIBS)




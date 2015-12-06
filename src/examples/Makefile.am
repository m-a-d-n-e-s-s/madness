include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt


noinst_PROGRAMS = pno h2dft hedft hello hatom_energy h2 he tdse heat heat2 csqrt \
 sdf_shape_tester test_gmres tdse1d vnucso nonlinschro sininteg functionio \
 dataloadbal hatom_1d binaryop dielectric hehf 3dharmonic testsolver \
 testspectralprop dielectric_external_field mp2 tiny oep h2dynamic newsolver \
 tdhf cc2 nemo helium_exact density_smoothing
 
if MADNESS_HAS_LIBXC
noinst_PROGRAMS += hefxc
endif

thisincludedir = $(includedir)/examples
thisinclude_HEADERS = molecularmask.h  nonlinsol.h spectralprop.h

AUTOMAKE_OPTIONS = subdir-objects

pno_SOURCES = pno.cc 

mp2_SOURCES = mp2.cc 

helium_exact_SOURCES = helium_exact.cc 

oep_SOURCES = oep.cc 

tiny_SOURCES = tiny.cc

dataloadbal_SOURCES = dataloadbal.cc

dielectric_external_field_SOURCES = dielectric_external_field.cc

testspectralprop_SOURCES = testspectralprop.cc

testsolver_SOURCES = testsolver.cc

tdhf_SOURCES = tdhf.cc

cc2_SOURCES = cc2.cc

nemo_SOURCES = nemo.cc 

density_smoothing_SOURCES = density_smoothing.cc

newsolver_SOURCES = newsolver.cc

hedft_SOURCES = hedft.cc 
hedft_LDADD = $(LIBCHEM) $(MRALIBS)

h2dft_SOURCES = h2dft.cc
h2dft_LDADD = $(LIBCHEM) $(MRALIBS)

if MADNESS_HAS_LIBXC
hefxc_SOURCES = hefxc.cc 
hefxc_LDADD = $(LIBCHEM) $(MRALIBS)
endif

sininteg_SOURCES = sininteg.cc

3dharmonic_SOURCES = 3dharmonic.cc

functionio_SOURCES = functionio.cc

nonlinschro_SOURCES = nonlinschro.cc

heat_SOURCES = heat.cc

heat2_SOURCES = heat2.cc

hello_SOURCES = hello.cc

dielectric_SOURCES = dielectric.cc

hatom_energy_SOURCES = hatom_energy.cc

h2_SOURCES = h2.cc
 
h2dynamic_SOURCES = h2dynamic.cc

he_SOURCES = he.cc

hehf_SOURCES = hehf.cc

# h2hf_SOURCES = h2hf.cc

tdse_SOURCES = tdse.cc

tdse1d_SOURCES = tdse1d.cc

csqrt_SOURCES = csqrt.cc

sdf_shape_tester_SOURCES = sdf_shape_tester.cc

test_gmres_SOURCES = test_gmres.cc

vnucso_SOURCES = vnucso.cc

hatom_1d_SOURCES = hatom_1d.cc

binaryop_SOURCES = binaryop.cc

LDADD = $(LIBCHEM) $(MRALIBS)


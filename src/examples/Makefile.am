include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

# specify the git commit version in the info functions
GITREV="MADNESS_GITREVISION=\"`git --git-dir=$(top_srcdir)/.git rev-parse HEAD`\" "

# default location for basis sets etc
AM_CPPFLAGS += -D$(GITREV)


noinst_PROGRAMS = pno h2dft hedft hello hatom_energy h2 he tdse_example heat heat2 csqrt \
 sdf_shape_tester test_gmres tdse1d vnucso nonlinschro sininteg functionio \
 dataloadbal hatom_1d binaryop dielectric hehf 3dharmonic testsolver \
 testspectralprop dielectric_external_field mp2 tiny oep h2dynamic newsolver \
 nemo helium_exact density_smoothing siam_example gaussian lrccs cc2
 
if MADNESS_HAS_LIBXC
noinst_PROGRAMS += hefxc
endif

thisincludedir = $(includedir)/examples
thisinclude_HEADERS = molecularmask.h  nonlinsol.h spectralprop.h

AUTOMAKE_OPTIONS = subdir-objects

pno_SOURCES = pno.cc 
pno_LDADD = $(LIBCHEM) $(MRALIBS)

siam_example_SOURCES = siam_example.cc 

mp2_SOURCES = mp2.cc 
mp2_LDADD = $(LIBCHEM) $(MRALIBS)

helium_exact_SOURCES = helium_exact.cc 
helium_exact_LDADD = $(LIBCHEM) $(MRALIBS)

oep_SOURCES = oep.cc 
oep_LDADD = $(LIBCHEM) $(MRALIBS)

tiny_SOURCES = tiny.cc
tiny_LDADD = $(LIBCHEM) $(MRALIBS)

dataloadbal_SOURCES = dataloadbal.cc

dielectric_external_field_SOURCES = dielectric_external_field.cc

testspectralprop_SOURCES = testspectralprop.cc

testsolver_SOURCES = testsolver.cc

lrccs_SOURCES = lrccs.cc
lrccs_LDADD = $(LIBCHEM) $(MRALIBS)

cc2_SOURCES = cc2.cc
cc2_LDADD = $(LIBCHEM) $(MRALIBS)
 
nemo_SOURCES = nemo.cc 
nemo_LDADD = $(LIBCHEM) $(MRALIBS)

density_smoothing_SOURCES = density_smoothing.cc
density_smoothing_LDADD = $(LIBCHEM) $(MRALIBS)

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

tdse_example_SOURCES = tdse_example.cc

tdse1d_SOURCES = tdse1d.cc

csqrt_SOURCES = csqrt.cc

sdf_shape_tester_SOURCES = sdf_shape_tester.cc

test_gmres_SOURCES = test_gmres.cc

vnucso_SOURCES = vnucso.cc

hatom_1d_SOURCES = hatom_1d.cc

binaryop_SOURCES = binaryop.cc

gaussian_SOURCES = gaussian.cc

LDADD = $(MRALIBS)
#LDADD = $(LIBCHEM) $(MRALIBS)


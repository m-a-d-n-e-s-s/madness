include $(top_srcdir)/config/MakeGlobal.am

SUBDIRS = tr1

bin_PROGRAMS = world testhashthreaded testdc objexample objexample2 testprof \
               helloworld taskexample testatomicint testhashdc array_test \
               gtest_example qtest worldprofile_test worldptr_test worldref_test \
               worldreduce_test functional_test

lib_LIBRARIES = libMADworld.a

thisincludedir = $(includedir)/world
thisinclude_HEADERS = archive.h print.h worldam.h worldfut.h worldmpi.h \
	worldtask.h array.h sharedptr.h worldgop.h worldobj.h bufar.h nodefaults.h \
	typestuff.h enable_if.h worlddep.h worldhash.h worldref.h worldtypes.h \
	dqueue.h parar.h vecar.h worldexc.h worldmem.h worldser.h worldthread.h \
	worldrmi.h safempi.h worldpapi.h worldmutex.h print_seq.h worldhashmap.h \
	worldrange.h atomicint.h posixmem.h worldptr.h deferred_cleanup.h world.h \
	worldprofile.h worldtime.h binfsar.h mpiar.h textfsar.h worlddc.h \
	deferred_deleter.h scopedptr.h worldreduce.h worldtask.h taskfn.h make_task.h \
	functional.h GPU_streams.h 

array_test_SOURCES = array_test.cc
array_test_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
array_test_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
array_test_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a #$(NVCCFLAGS)

gtest_example_SOURCES = gtest_example.cc
gtest_example_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
gtest_example_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
gtest_example_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a $(NVCCFLAGS)

worldptr_test_SOURCES = worldptr_test.cc
worldptr_test_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
worldptr_test_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
worldptr_test_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a $(NVCCFLAGS)

worldref_test_SOURCES = worldref_test.cc
worldref_test_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
worldref_test_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
worldref_test_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a  $(NVCCFLAGS)

worldreduce_test_SOURCES = worldreduce_test.cc
worldreduce_test_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
worldreduce_test_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
worldreduce_test_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a $(NVCCFLAGS)

functional_test_SOURCES = functional_test.cc
functional_test_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
functional_test_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
functional_test_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a $(NVCCFLAGS)

testprof_SOURCES = testprof.cc
testprof_LDADD = libMADworld.a $(NVCCFLAGS)

testhashdc_SOURCES = testhashdc.cc
testhashdc_LDADD = libMADworld.a $(NVCCFLAGS)

taskexample_SOURCES = taskexample.cc
taskexample_LDADD = libMADworld.a $(NVCCFLAGS)

helloworld_SOURCES = helloworld.cc
helloworld_LDADD = libMADworld.a $(NVCCFLAGS)

testatomicint_SOURCES = testatomicint.cc
testatomicint_LDADD = libMADworld.a $(NVCCFLAGS)

objexample_SOURCES = objexample.cc
objexample_LDADD = libMADworld.a $(NVCCFLAGS)

objexample2_SOURCES = objexample2.cc
objexample2_LDADD = libMADworld.a $(NVCCFLAGS)

testdc_SOURCES = testdc.cc
testdc_LDADD = libMADworld.a $(NVCCFLAGS)

testhashthreaded_SOURCES = testhashthreaded.cc
testhashthreaded_LDADD = libMADworld.a  $(NVCCFLAGS)

qtest_SOURCES = qtest.cc 
qtest_LDADD = libMADworld.a $(NVCCFLAGS)

if HAVE_CUDA
world_SOURCES = world.cc cuda_streams.cu
else
world_SOURCES = world.cc
endif

world_LDADD = libMADworld.a $(NVCCFLAGS)

worldprofile_test_SOURCES = worldprofile_test.cc
worldprofile_test_LDADD = libMADworld.a $(NVCCFLAGS)

if HAVE_CUDA
libMADworld_a_SOURCES =cuda_streams.cu  worldstuff.cc redirectio.cc archive_type_names.cc \
	debug.cc print.cc worldmem.cc lookup3.c worldrmi.cc safempi.cc worldpapi.cc \
	worldref.cc worldam.cc worldprofile.cc worldthread.cc worldtask.cc \
	worldgop.cc deferred_cleanup.cc worlddep.cc worldmutex.cc binfsar.cc \
	worldmem.cc textfsar.cc \
	$(thisinclude_HEADERS)

.cu.o: cuda_stream.cu
	$(NVCC)  -arch=sm_20 --compiler-options="-DHAVE_CONFIG_H   -I../../../include -I../../../src/lib  -I../../../src/apps -I../../../include  -g -Wall -Wno-strict-aliasing -Wno-deprecated  -ffast-math -march=native -mfpmath=sse -msse -mpc64" $(NVCCFLAGS) $(AM_CPPFLAGS) -c  -o $@ $<
libraries:	$(lib_LIBRARIES)
else
libMADworld_a_SOURCES =worldstuff.cc redirectio.cc archive_type_names.cc \
	debug.cc print.cc worldmem.cc lookup3.c worldrmi.cc safempi.cc worldpapi.cc \
	worldref.cc worldam.cc worldprofile.cc worldthread.cc worldtask.cc \
	worldgop.cc deferred_cleanup.cc worlddep.cc worldmutex.cc binfsar.cc \
	worldmem.cc textfsar.cc \
	$(thisinclude_HEADERS)
endif

	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

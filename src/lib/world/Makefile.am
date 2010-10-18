include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = world testhashthreaded testdc objexample objexample2 testprof helloworld taskexample testatomicint testhashdc array_test world_test qtest worldprofile_test

lib_LIBRARIES = libMADworld.a

thisincludedir = $(includedir)/world
thisinclude_HEADERS = archive.h print.h      worldam.h       worldfut.h   worldmpi.h        worldtask.h \
                     array.h    sharedptr.h  worldgop.h   worldobj.h \
                     binfsar.h  mpiar.h       textfsar.h   worlddc.h       world.h      worldprofile.h  worldtime.h \
                     bufar.h    nodefaults.h  typestuff.h  worlddep.h      worldhash.h  worldref.h      worldtypes.h \
                     dqueue.h   parar.h       vecar.h      worldexc.h      worldmem.h   worldser.h \
                     worldthread.h worldrmi.h safempi.h worldpapi.h worldmutex.h print_seq.h \
                     worldhashmap.h          worldrange.h  atomicint.h     posixmem.h

array_test_SOURCES = array_test.cc
array_test_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
array_test_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
array_test_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS)

world_test_SOURCES = world_test.cc
world_test_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
world_test_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
world_test_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a

testprof_SOURCES = testprof.cc
testprof_LDADD = libMADworld.a

testhashdc_SOURCES = testhashdc.cc
testhashdc_LDADD = libMADworld.a

taskexample_SOURCES = taskexample.cc
taskexample_LDADD = libMADworld.a

helloworld_SOURCES = helloworld.cc
helloworld_LDADD = libMADworld.a

testatomicint_SOURCES = testatomicint.cc threadpool.cc
testatomicint_LDADD = libMADworld.a

objexample_SOURCES = objexample.cc threadpool.cc
objexample_LDADD = libMADworld.a

objexample2_SOURCES = objexample2.cc threadpool.cc
objexample2_LDADD = libMADworld.a

testdc_SOURCES = testdc.cc
testdc_LDADD = libMADworld.a

testhashthreaded_SOURCES = testhashthreaded.cc threadpool.cc
testhashthreaded_LDADD = libMADworld.a

qtest_SOURCES = qtest.cc
qtest_LDADD = libMADworld.a


world_SOURCES = world.cc
world_LDADD = libMADworld.a

worldprofile_test_SOURCES = worldprofile_test.cc
worldprofile_test_LDADD = libMADworld.a

libMADworld_a_SOURCES = worldstuff.cc redirectio.cc archive_type_names.cc debug.cc print.cc worldmem.cc lookup3.c \
                         threadpool.cc worldrmi.cc safempi.cc worldpapi.cc worldref.cc \
                         worldam.cc worldprofile.cc worldthread.cc worldtask.cc worldgop.cc deferred_cleanup.cc \
                         worlddep.cc worldmutex.cc \
                     \
                     $(thisinclude_HEADERS)


libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

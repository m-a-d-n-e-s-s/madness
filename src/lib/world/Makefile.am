include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test_world test_hashthreaded test_dc test_future2 test_future3 test_prof \
               test_hello test_future test_atomicint test_hashdc test_array test_ar \
               test_googletest test_queue test_worldprofile test_worldptr test_worldref

lib_LIBRARIES = libMADworld.a

thisincludedir = $(includedir)/world
thisinclude_HEADERS = archive.h print.h worldam.h worldfut.h worldmpi.h \
	worldtask.h array.h worldgop.h worldobj.h bufar.h nodefaults.h \
	typestuff.h enable_if.h worlddep.h worldhash.h worldref.h worldtypes.h \
	dqueue.h parar.h vecar.h worldexc.h worldmem.h worldser.h worldthread.h \
	worldrmi.h safempi.h worldpapi.h worldmutex.h print_seq.h worldhashmap.h \
	worldrange.h atomicint.h posixmem.h worldptr.h deferred_cleanup.h world.h \
	worldfwd.h uniqueid.h worldprofile.h worldtime.h binfsar.h mpiar.h textfsar.h \
	worlddc.h deferred_deleter.h scopedptr.h worldtask.h taskfn.h \
	shared_ptr.h shared_ptr_bits.h stdarray.h ref.h move.h group.h dist_cache.h dist_keys.h \
	stdarray_bits.h type_traits.h type_traits_bits.h boost_checked_delete_bits.h \
	function_traits_bits.h integral_constant.h stubmpi.h bgq_atomics.h
                      
TESTS = test_array test_worldptr test_worldref test_prof \
        test_hashdc test_hello test_atomicint test_future2 \
        test_future3 test_dc test_hashthreaded test_queue test_world test_worldprofile 
                      
test_array_SOURCES = test_array.cc
test_array_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
test_array_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
test_array_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a

test_googletest_SOURCES = test_googletest.cc
test_googletest_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
test_googletest_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
test_googletest_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a

test_worldptr_SOURCES = test_worldptr.cc
test_worldptr_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
test_worldptr_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
test_worldptr_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a

test_worldref_SOURCES = test_worldref.cc
test_worldref_CPPFLAGS = $(GTEST_CPPFLAGS) $(AM_CPPFLAGS)
test_worldref_CXXFLAGS = $(GTEST_CXXFLAGS) $(AM_CXXFLAGS)
test_worldref_LDADD = $(GTEST_LDFLAGS) $(GTEST_LIBS) libMADworld.a

test_prof_SOURCES = test_prof.cc
test_prof_LDADD = libMADworld.a

test_ar_SOURCES = test_ar.cc
test_ar_LDADD = libMADworld.a

test_hashdc_SOURCES = test_hashdc.cc
test_hashdc_LDADD = libMADworld.a

test_future_SOURCES = test_future.cc
test_future_LDADD = libMADworld.a

test_hello_SOURCES = test_hello.cc
test_hello_LDADD = libMADworld.a

test_atomicint_SOURCES = test_atomicint.cc
test_atomicint_LDADD = libMADworld.a

test_future2_SOURCES = test_future2.cc
test_future2_LDADD = libMADworld.a

test_future3_SOURCES = test_future3.cc
test_future3_LDADD = libMADworld.a

test_dc_SOURCES = test_dc.cc
test_dc_LDADD = libMADworld.a

test_hashthreaded_SOURCES = test_hashthreaded.cc
test_hashthreaded_LDADD = libMADworld.a

test_queue_SOURCES = test_queue.cc
test_queue_LDADD = libMADworld.a

test_world_SOURCES = test_world.cc
test_world_LDADD = libMADworld.a

test_worldprofile_SOURCES = test_worldprofile.cc
test_worldprofile_LDADD = libMADworld.a

libMADworld_a_SOURCES = world.cc redirectio.cc archive_type_names.cc \
	debug.cc print.cc worldmem.cc worldrmi.cc safempi.cc worldpapi.cc \
	worldref.cc worldam.cc worldprofile.cc worldthread.cc worldtask.cc \
	worldgop.cc deferred_cleanup.cc worldmutex.cc binfsar.cc textfsar.cc \
    lookup3.c worldmpi.cc group.cc \
	$(thisinclude_HEADERS)


libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

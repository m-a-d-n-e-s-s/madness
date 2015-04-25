include $(top_srcdir)/config/MakeGlobal.am


lib_LIBRARIES = libMADworld.a

thisincludedir = $(includedir)/madness/world
thisinclude_HEADERS = archive.h print.h worldam.h future.h worldmpi.h \
	world_task_queue.h array.h worldgop.h world_object.h buffer_archive.h \
	nodefaults.h enable_if.h worlddep.h worldhash.h worldref.h worldtypes.h \
	dqueue.h parar.h vector_archive.h madness_exception.h worldmem.h \
	thread.h worldrmi.h safempi.h worldpapi.h worldmutex.h print_seq.h \
	worldhashmap.h worldrange.h atomicint.h posixmem.h worldptr.h \
	deferred_cleanup.h parallel_runtime.h world.h uniqueid.h worldprofile.h \
	timers.h binary_fstream_archive.h mpiar.h text_fstream_archive.h worlddc.h \
	mem_func_wrapper.h scopedptr.h taskfn.h group.h dist_cache.h \
	dist_keys.h type_traits.h boost_checked_delete_bits.h \
	function_traits.h stubmpi.h bgq_atomics.h binsorter.h


                      
TESTS = test_prof.mpi test_ar.mpi test_hashdc.mpi test_hello.mpi test_atomicint.mpi test_future.mpi \
        test_future2.mpi test_future3.mpi test_dc.mpi test_hashthreaded.mpi test_queue.mpi test_world.mpi \
        test_worldprofile.mpi test_binsorter.mpi


if MADNESS_HAS_GOOGLE_TEST
TESTS += test_array.mpi test_worldptr.mpi test_worldref.mpi
XFAIL_TESTS =  test_googletest.mpi
endif

TEST_EXTENSIONS = .mpi .seq

# tests run by mpirun
MPI_LOG_COMPILER = $(MPIRUN)
AM_MPI_LOG_FLAGS = $(MPIFLAGS)

# sequential tests
SEQ_LOG_COMPILER = 
AM_SEQ_LOG_FLAGS =

# default test (= sequential)
LOG_COMPILER = 
AM_LOG_FLAGS =




noinst_PROGRAMS = $(TESTS)

test_prof_mpi_SOURCES = test_prof.cc
test_prof_mpi_LDADD = libMADworld.a

test_binsorter_mpi_SOURCES = test_binsorter.cc
test_binsorter_mpi_LDADD = libMADworld.a

test_ar_mpi_SOURCES = test_ar.cc
test_ar_mpi_LDADD = libMADworld.a

test_hashdc_mpi_SOURCES = test_hashdc.cc
test_hashdc_mpi_LDADD = libMADworld.a

test_future_mpi_SOURCES = test_future.cc
test_future_mpi_LDADD = libMADworld.a

test_hello_mpi_SOURCES = test_hello.cc
test_hello_mpi_LDADD = libMADworld.a

test_atomicint_mpi_SOURCES = test_atomicint.cc
test_atomicint_mpi_LDADD = libMADworld.a

test_future2_mpi_SOURCES = test_future2.cc
test_future2_mpi_LDADD = libMADworld.a

test_future3_mpi_SOURCES = test_future3.cc
test_future3_mpi_LDADD = libMADworld.a

test_dc_mpi_SOURCES = test_dc.cc
test_dc_mpi_LDADD = libMADworld.a

test_hashthreaded_mpi_SOURCES = test_hashthreaded.cc
test_hashthreaded_mpi_LDADD = libMADworld.a

test_queue_mpi_SOURCES = test_queue.cc
test_queue_mpi_LDADD = libMADworld.a

test_world_mpi_SOURCES = test_world.cc
test_world_mpi_LDADD = libMADworld.a

test_worldprofile_mpi_SOURCES = test_worldprofile.cc
test_worldprofile_mpi_LDADD = libMADworld.a

if MADNESS_HAS_GOOGLE_TEST

test_array_mpi_SOURCES = test_array.cc
test_array_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_array_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_array_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.a

#test_googletest_mpi_SOURCES = test_googletest.cc
#test_googletest_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
#test_googletest_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
#test_googletest_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.a

test_worldptr_mpi_SOURCES = test_worldptr.cc
test_worldptr_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_worldptr_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_worldptr_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.a

test_worldref_mpi_SOURCES = test_worldref.cc
test_worldref_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_worldref_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_worldref_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.a

endif

libMADworld_a_SOURCES = madness_exception.cc world.cc timers.cc future.cc \
	redirectio.cc archive_type_names.cc \
	debug.cc print.cc worldmem.cc worldrmi.cc safempi.cc worldpapi.cc \
	worldref.cc worldam.cc worldprofile.cc thread.cc world_task_queue.cc \
	worldgop.cc deferred_cleanup.cc worldmutex.cc binary_fstream_archive.cc \
	text_fstream_archive.cc lookup3.c worldmpi.cc group.cc \
	$(thisinclude_HEADERS)


libraries:	$(lib_LIBRARIES)

install-libraries: install-libLIBRARIES install-thisincludeHEADERS

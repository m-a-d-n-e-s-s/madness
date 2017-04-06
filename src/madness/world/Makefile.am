include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

# specify the git commit version in the info functions
GITREV="MADNESS_GITREVISION=\"`git --git-dir=$(top_srcdir)/.git rev-parse HEAD`\" "

lib_LTLIBRARIES = libMADworld.la

thisincludedir = $(includedir)/madness/world
thisinclude_HEADERS = info.h archive.h print.h worldam.h future.h worldmpi.h \
	world_task_queue.h array_addons.h stack.h vector.h worldgop.h \
	world_object.h buffer_archive.h \
	nodefaults.h dependency_interface.h worldhash.h worldref.h worldtypes.h \
	dqueue.h parallel_archive.h vector_archive.h madness_exception.h \
	worldmem.h thread.h worldrmi.h safempi.h worldpapi.h worldmutex.h \
	print_seq.h worldhashmap.h range.h atomicint.h posixmem.h worldptr.h \
	deferred_cleanup.h MADworld.h world.h uniqueid.h worldprofile.h \
	timers.h binary_fstream_archive.h mpi_archive.h text_fstream_archive.h \
	worlddc.h mem_func_wrapper.h taskfn.h group.h dist_cache.h \
	distributed_id.h type_traits.h \
	function_traits.h stubmpi.h bgq_atomics.h binsorter.h


                      
TESTS = test_prof.mpi test_ar.mpi test_hashdc.mpi test_hello.mpi test_atomicint.mpi test_future.mpi \
        test_future2.mpi test_future3.mpi test_dc.mpi test_hashthreaded.mpi test_queue.mpi test_world.mpi \
        test_worldprofile.mpi test_binsorter.mpi test_tree.mpi


if MADNESS_HAS_GOOGLE_TEST
TESTS += test_vector.mpi test_worldptr.mpi test_worldref.mpi test_stack.seq
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
test_prof_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_binsorter_mpi_SOURCES = test_binsorter.cc
test_binsorter_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_ar_mpi_SOURCES = test_ar.cc
test_ar_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_hashdc_mpi_SOURCES = test_hashdc.cc
test_hashdc_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_future_mpi_SOURCES = test_future.cc
test_future_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_hello_mpi_SOURCES = test_hello.cc
test_hello_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_atomicint_mpi_SOURCES = test_atomicint.cc
test_atomicint_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_future2_mpi_SOURCES = test_future2.cc
test_future2_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_future3_mpi_SOURCES = test_future3.cc
test_future3_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_dc_mpi_SOURCES = test_dc.cc
test_dc_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_tree_mpi_SOURCES = test_tree.cc
test_tree_mpi_LDADD = libMADworld.la

test_hashthreaded_mpi_SOURCES = test_hashthreaded.cc
test_hashthreaded_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_queue_mpi_SOURCES = test_queue.cc
test_queue_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_world_mpi_SOURCES = test_world.cc
test_world_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

test_worldprofile_mpi_SOURCES = test_worldprofile.cc
test_worldprofile_mpi_LDADD = libMADworld.la ${PaRSEC_LIBS}

if MADNESS_HAS_GOOGLE_TEST

test_vector_mpi_SOURCES = test_vector.cc
test_vector_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_vector_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_vector_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.la

#test_googletest_mpi_SOURCES = test_googletest.cc
#test_googletest_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
#test_googletest_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
#test_googletest_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.la

test_worldptr_mpi_SOURCES = test_worldptr.cc
test_worldptr_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_worldptr_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_worldptr_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.la

test_worldref_mpi_SOURCES = test_worldref.cc
test_worldref_mpi_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_worldref_mpi_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_worldref_mpi_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.la

test_stack_seq_SOURCES = test_stack.cc
test_stack_seq_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_stack_seq_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_stack_seq_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADworld.la

endif

libMADworld_la_SOURCES = madness_exception.cc world.cc timers.cc future.cc \
	redirectio.cc archive_type_names.cc info.cc \
	debug.cc print.cc worldmem.cc worldrmi.cc safempi.cc worldpapi.cc \
	worldref.cc worldam.cc worldprofile.cc thread.cc world_task_queue.cc \
	worldgop.cc deferred_cleanup.cc worldmutex.cc binary_fstream_archive.cc \
	text_fstream_archive.cc lookup3.c worldmpi.cc group.cc \
	$(thisinclude_HEADERS)

libMADworld_la_CPPFLAGS = $(AM_CPPFLAGS) -D$(GITREV)
libMADworld_la_LDFLAGS = -version-info 0:0:0
if USE_PARSEC
libMADworld_la_SOURCES += parsec.cc
endif


libraries:	$(lib_LTLIBRARIES)

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS

MOSTLYCLEANFILES = test.dat testme.ar

include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = world testatom testhashthreaded testdc objexample objexample2 testprof helloworld
lib_LIBRARIES = libMADworld.a

thisincludedir = $(includedir)/world
thisinclude_HEADERS = archive.h print.h      worldam.h       worldfut.h   worldmpi.h        worldtask.h \
                     array.h    sharedptr.h  worldgop.h   worldobj.h \
                     binfsar.h  mpiar.h       textfsar.h   worlddc.h       world.h      worldprofile.h  worldtime.h \
                     bufar.h    nodefaults.h  typestuff.h  worlddep.h      worldhash.h  worldref.h      worldtypes.h \
                     dqueue.h   parar.h       vecar.h      worldexc.h      worldmem.h   worldser.h \
                     worldthread.h madatomic.h worldrmi.h safempi.h worldpapi.h worldmutex.h print_seq.h \
                     worldhashmap.h          worldrange.h

testprof_SOURCES = testprof.cc
testprof_LDADD = libMADworld.a

helloworld_SOURCES = helloworld.cc
helloworld_LDADD = libMADworld.a

testatom_SOURCES = testatom.cc threadpool.cc
testatom_LDADD = libMADworld.a

objexample_SOURCES = objexample.cc threadpool.cc
objexample_LDADD = libMADworld.a

objexample2_SOURCES = objexample2.cc threadpool.cc
objexample2_LDADD = libMADworld.a

testdc_SOURCES = testdc.cc
testdc_LDADD = libMADworld.a

testhashthreaded_SOURCES = testhashthreaded.cc threadpool.cc
testhashthreaded_LDADD = libMADworld.a


world_SOURCES = world.cc
world_LDADD = libMADworld.a

libMADworld_a_SOURCES = worldstuff.cc redirectio.cc archive_type_names.cc debug.cc print.cc worldmem.cc lookup3.c \
                         threadpool.cc worldrmi.cc safempi.cc worldpapi.cc \
                     \
                     archive.h print.h      worldam.h       worldfut.h   worldmpi.h        worldtask.h \
                     array.h    sharedptr.h  worldgop.h   worldobj.h \
                     binfsar.h  mpiar.h       textfsar.h   worlddc.h       world.h      worldprofile.h  worldtime.h \
                     bufar.h    nodefaults.h  typestuff.h  worlddep.h      worldhash.h  worldref.h      worldtypes.h \
                     dqueue.h   parar.h       vecar.h      worldexc.h      worldmem.h   worldser.h \
                     worldthread.h madatomic.h worldrmi.h safempi.h worldpapi.h worldmutex.h print_seq.h



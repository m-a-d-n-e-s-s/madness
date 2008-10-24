include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = world testatom testhashthreaded
lib_LIBRARIES = libMADworld.a

testatom_SOURCES = testatom.cc threadpool.cc

testhashthreaded_SOURCES = testhashthreaded.cc threadpool.cc
testhashthreaded_LDADD = libMADworld.a


world_SOURCES = world.cc
world_LDADD = libMADworld.a

libMADworld_a_SOURCES = worldstuff.cc redirectio.cc archive_type_names.cc debug.cc print.cc worldmem.cc lookup3.c \
                         threadpool.cc worldrmi.cc safempi.cc \
                     \
                     archive.h print.h      worldam.h       worldfut.h   worldmpi.h        worldtask.h \
                     array.h    sharedptr.h  worldgop.h   worldobj.h \
                     binfsar.h  mpiar.h       textfsar.h   worlddc.h       world.h      worldprofile.h  worldtime.h \
                     bufar.h    nodefaults.h  typestuff.h  worlddep.h      worldhash.h  worldref.h      worldtypes.h \
                     dqueue.h   parar.h       vecar.h      worldexc.h      worldmem.h   worldser.h \
                     worldthread.h madatomic.h worldrmi.h safempi.h


include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS    = tdselb bigboy2lb bigboylb testshqi navstokes_cosines_chorin navstokes_cosines_lb_iter 

tdselb_SOURCES = tdselb.cc
tdselb_LDADD   = $(MRALIBS)

bigboy2lb_SOURCES = bigboy2lb.cc
bigboy2lb_LDADD   = $(MRALIBS)

bigboylb_SOURCES = bigboylb.cc molecule.h molecule.cc
bigboylb_LDADD   = $(MRALIBS)

testshqi_SOURCES = testshqi.cc
testshqi_LDADD   = $(MRALIBS)

navstokes_cosines_chorin_SOURCES = navstokes_cosines_chorin.cc poperator.h mybicgstab.h
navstokes_cosines_chorin_LDADD = $(MRALIBS)

navstokes_cosines_lb_iter_SOURCES = navstokes_cosines_lb_iter.cc poperator.h mybicgstab.h solutions.h
navstokes_cosines_lb_iter_LDADD = $(MRALIBS)

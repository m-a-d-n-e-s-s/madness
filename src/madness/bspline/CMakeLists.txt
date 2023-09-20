set (MADBSPLINE_SOURCES
  bspline.cc
  )

set (MADBSPLINE_HEADERS
  gauleg.h
  myqd.h
  )

add_mad_executable(bspline bspline.cc "MADchem;-lqd")
##add_dependencies(applications-madness bspline)



set (MADBSPLINE_SOURCES
  bspline.cc madbessel.cc solidharmonics.cc
  )

set (MADBSPLINE_HEADERS
  gauleg.h
  myqd.h
  )

add_mad_executable(bspline bspline.cc "MADchem;-lqd")
add_mad_executable(solidharmonics solidharmonics.cc "MADchem;-lqd")
add_mad_executable(madbessel madbessel.cc "MADchem;-lqd")
##add_dependencies(applications-madness bspline)



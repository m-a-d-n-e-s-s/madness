
#Find HDF5 support
if (ENABLE_HDF5)
  set(HDF5_PREFER_PARALLEL FALSE)
  find_package(HDF5 COMPONENTS C)

  if (HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_LIBRARIES}")
    set(HAVE_HDF5 ON CACHE INTERNAL "Have HDF5" FORCE)
  endif ()
#  if (HDF5_IS_PARALLEL)
#    set(HAVE_HDF5 0)
#    message(STATUS "Found HDF5 Parallel -  Disable")
#    list(REMOVE_ITEM HDF5_C_LIBRARY_NAMES)
#  endif ()

endif ()

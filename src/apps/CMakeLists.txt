# src/apps

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  
add_subdirectory(chem)
add_subdirectory(tdse)
add_subdirectory(moldft)
add_subdirectory(molresponse)
add_subdirectory(pno)
add_subdirectory(cis)
#add_subdirectory(hf)
#add_subdirectory(interior_bc)
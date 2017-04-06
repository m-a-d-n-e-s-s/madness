# doc

if(DOXYGEN_FOUND)
  set_directory_properties(PROPERTY ADDITIONAL_MAKE_CLEAN_FILES html latex *~)

  add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} doxygen.cfg)
  
  install(DIRECTORY html 
      DESTINATION "${MADNESS_INSTALL_DOCDIR}" OPTIONAL
      COMPONENT docs)
endif()
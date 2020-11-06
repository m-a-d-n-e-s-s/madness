# doc

if(DOXYGEN_FOUND)
  set_directory_properties(PROPERTY ADDITIONAL_MAKE_CLEAN_FILES html latex *~)

  add_custom_target_subproject(madness html
    COMMAND ${DOXYGEN_EXECUTABLE} doxygen.cfg
    USES_TERMINAL)

  add_custom_target_subproject(madness docs)
  add_dependencies(docs-madness html-madness)

  add_dependencies(everything docs)
  add_dependencies(everything docs-madness)
  add_dependencies(everything html-madness)

  install(DIRECTORY html
      DESTINATION "${MADNESS_INSTALL_DOCDIR}" OPTIONAL
      COMPONENT docs)
endif()

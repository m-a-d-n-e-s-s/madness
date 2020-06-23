if(ENABLE_PARSEC)
  if (NOT TARGET PaRSEC::parsec)
    find_package(PARSEC REQUIRED)
  endif(NOT TARGET PaRSEC::parsec)

  if(TARGET PaRSEC::parsec)
    set(HAVE_PARSEC 1)
  else()
    message(FATAL_ERROR "PaRSEC found but PaRSEC::parsec target not provided; using legacy PaRSEC build?")
  endif()
endif()

# if target PaRSEC::parsec is an ALIAS to parsec, we are using parsec as a subproject, need to add parsec to export sets for all dependent targets
macro(add_subproject_parsec_target_to_export_set _export_set _component)

  if (TARGET PaRSEC::parsec AND TARGET parsec)
    get_property(_parsecparsec_aliased_target_is_set TARGET PaRSEC::parsec PROPERTY ALIASED_TARGET SET)
    if (_parsecparsec_aliased_target_is_set)  # PaRSEC::parsec is an ALIAS target ... probably to parsec, but check
      get_property(_parsecparsec_aliased_target TARGET PaRSEC::parsec PROPERTY ALIASED_TARGET)
      if (NOT (_parsecparsec_aliased_target STREQUAL parsec) )
        message(FATAL_ERROR "TARGET \"PaRSEC::parsec\" is an alias, but is aliased to target other than \"parsec\"")
      endif()

      install(TARGETS parsec EXPORT "${_export_set}"
          COMPONENT "${_component}"
          LIBRARY DESTINATION "${MADNESS_INSTALL_LIBDIR}"
          ARCHIVE DESTINATION "${MADNESS_INSTALL_LIBDIR}"
          INCLUDES DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}")

    endif()
  endif()

endmacro()

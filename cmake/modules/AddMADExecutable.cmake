macro(add_mad_executable _name _source_files _libs)

  add_executable(${_name} ${_source_files})
  target_link_libraries(${_name} PRIVATE "${_libs}")
  add_dependencies(everything ${_name})

endmacro()

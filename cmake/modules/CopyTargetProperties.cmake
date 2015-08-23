macro(copy_target_properties _intarget _outtarget _props)

  foreach(_prop ${_props})
    get_property(_prop_set TARGET ${_intarget} PROPERTY ${_prop} SET)
    if(_prop_set)
      get_target_property(_value ${_intarget} ${_prop})
      set_target_properties(${_outtarget} PROPERTIES ${_prop} "${_value}")
    endif()
  endforeach()

endmacro()
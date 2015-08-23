macro(copy_target_properties _intarget _outtarget _props)

  foreach(_prop ${_props})
    get_target_property(_value ${_intarget} ${_prop})
    if(_value)
      set_target_properties(${_outtarget} PROPERTIES ${_prop} "${_value}")
    endif()
  endforeach()

endmacro()
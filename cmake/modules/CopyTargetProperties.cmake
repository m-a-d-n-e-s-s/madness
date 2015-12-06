macro(copy_target_properties _intarget _outtarget _props)

  foreach(_prop ${_props})
    get_property(_prop_set TARGET ${_intarget} PROPERTY ${_prop} SET)
    get_property(_prop_value TARGET ${_intarget} PROPERTY ${_prop})
#    message("${_intarget}:${_prop} SET = '${_prop_set}'")
#    message("${_intarget}:${_prop} = '${_prop_value}'")
    if(_prop_set)
      set_property(TARGET ${_outtarget} PROPERTY ${_prop} "${_prop_value}")
      
#      get_target_property(_new_value ${_outtarget} ${_prop})
#      message("${_outtarget}:${_prop} = '${_new_value}'")
    endif()
    
    unset(_prop_set)
    unset(_prop_value)
  endforeach()

endmacro()

macro(append_target_properties _intarget _outtarget _props)

  foreach(_prop ${_props})
    get_property(_prop_set TARGET ${_intarget} PROPERTY ${_prop} SET)
    get_property(_prop_value TARGET ${_intarget} PROPERTY ${_prop})
#    message("${_intarget}:${_prop} SET = '${_prop_set}'")
#    message("${_intarget}:${_prop} = '${_prop_value}'")
    if(_prop_set)
      set_property(TARGET ${_outtarget} APPEND PROPERTY ${_prop} "${_prop_value}")
      
#      get_target_property(_new_value ${_outtarget} ${_prop})
#      message("${_outtarget}:${_prop} = '${_new_value}'")
    endif()
    
    unset(_prop_set)
    unset(_prop_value)
  endforeach()

endmacro()
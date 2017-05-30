
include(CheckCXXCompilerFlag)

macro(check_disablepie_support _outvar _disablepie_linker_flags)

  if(NOT ${_outvar})
    message(STATUS "Checking for PIE-disabling linker flags")
  endif()
  
  # set the flag manually for Darwin
  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(disablepie_linker_flags "-Wl,-no_pie")
  else()
    set(disablepie_linker_flags )
    foreach(_disablepie_test_flag "-no-pie")
      
      # Try compiling
      unset(${_outvar} CACHE)
      check_cxx_compiler_flag(${_disablepie_test_flag} ${_outvar})
      
      if(${_outvar})
        list(APPEND disablepie_linker_flags "${_disablepie_test_flag}")
        break()
      endif()
      
    endforeach()
  endif()

  if (disablepie_linker_flags)
    set(${_disablepie_linker_flags} "${disablepie_linker_flags}"
          CACHE STRING "Linker flags required to disable PIE support")
    mark_as_advanced(${_disablepie_linker_flags})
    message(STATUS "PIE-disabling linker flags: ${${_disablepie_linker_flags}}")
  endif()
  
endmacro(check_disablepie_support)
if (ENABLE_TBB)
  find_package(TBB 4.3.5)

  if (TBB_FOUND AND TARGET TBB::tbb)
    set(HAVE_INTEL_TBB 1)
    message(STATUS "Found TBB: TBB_CONFIG=${TBB_CONFIG}")
  endif ()

  # if user set MADNESS_EXPLOIT_TBB_PRIORITY to ON:
  # check if task priorities are supported (albeit deprecated) and whether they need macro __TBB_TASK_PRIORITY to be defined
  if (MADNESS_EXPLOIT_TBB_PRIORITY)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_INCLUDES ${TBB_INCLUDE_DIRS})
    if (TBB_USE_DEBUG AND TBB_LIBRARIES_DEBUG)
      list(APPEND CMAKE_REQUIRED_DEFINITIONS "${TBB_COMPILE_FLAGS_DEBUG}")
      set(CMAKE_REQUIRED_LIBRARIES ${TBB_LIBRARIES_DEBUG})
    else ()
      set(CMAKE_REQUIRED_LIBRARIES ${TBB_LIBRARIES_RELEASE})
    endif ()
    set(_tbb_task_priority_compile_check_source "
#include <tbb/task.h>
int main(int argc, char** argv) {
tbb::task* t;
tbb::task::enqueue(*t, tbb::priority_high);
}
")

    check_cxx_source_compiles("${_tbb_task_priority_compile_check_source}" MADNESS_CAN_USE_TBB_PRIORITY)
    if (MADNESS_CAN_USE_TBB_PRIORITY)
      set(MADNESS_CAN_USE_TBB_PRIORITY TRUE CACHE BOOL "Whether MADNESS can control TBB task priorities" FORCE)
    else (MADNESS_CAN_USE_TBB_PRIORITY)
      list(APPEND CMAKE_REQUIRED_DEFINITIONS "-D__TBB_TASK_PRIORITY=1")
      check_cxx_source_compiles("${_tbb_task_priority_compile_check_source}" MADNESS_CAN_USE_TBB_PRIORITY_IF_MACROED_ON)
      if (MADNESS_CAN_USE_TBB_PRIORITY_IF_MACROED_ON)
        set(MADNESS_CAN_USE_TBB_PRIORITY TRUE CACHE BOOL "Whether MADNESS can control TBB task priorities" FORCE)
        list(APPEND TBB_COMPILE_FLAGS "-D__TBB_TASK_PRIORITY=1")
      endif (MADNESS_CAN_USE_TBB_PRIORITY_IF_MACROED_ON)
    endif (MADNESS_CAN_USE_TBB_PRIORITY)
    if (MADNESS_CAN_USE_TBB_PRIORITY)
      message(STATUS "TBB task priorities will be enabled")
    endif (MADNESS_CAN_USE_TBB_PRIORITY)

    cmake_pop_check_state()
  endif (MADNESS_EXPLOIT_TBB_PRIORITY)


endif ()
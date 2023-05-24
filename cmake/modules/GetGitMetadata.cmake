#
# SPDX-FileCopyrightText: 2022 Eduard Valeyev <eduard@valeyev.net>
#
# SPDX-License-Identifier: BSD-2-Clause
#

#
# vgkit_cmake_git_revision() defines variable ${PROJECT_NAME_UPPER}_GIT_REVISION,
#   where ${PROJECT_NAME_UPPER} is the upper-case version of ${PROJECT_NAME}, to
#   the "long-form" GIT revision produced by `git rev-parse -q HEAD` command.
# N.B. If the variable is already defined will print a warning message but not change the value of
#   ${PROJECT_NAME_UPPER}_GIT_REVISION
#
macro(vgkit_cmake_git_revision)
  string(TOUPPER "${PROJECT_NAME}_GIT_REVISION" vgkit_cmake_git_revision_result)
  if (DEFINED ${vgkit_cmake_git_revision_result})
    message(WARNING "vgkit_cmake_git_revision() called for project ${PROJECT_NAME} but the result variable ${vgkit_cmake_git_revision_result} already defined; will not change its value")
  else()
    if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
      find_package(Git REQUIRED)
      if (GIT_FOUND)
        execute_process(
          COMMAND ${GIT_EXECUTABLE} rev-parse -q HEAD
          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
          RESULT_VARIABLE vgkit_cmake_git_revision_errcod
          OUTPUT_VARIABLE ${vgkit_cmake_git_revision_result}
          ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT vgkit_cmake_git_revision_errcod EQUAL 0)
          set(${vgkit_cmake_git_revision_result} "unknown")
        endif()
        unset(vgkit_cmake_git_revision_errcod)
      else()
        set(${vgkit_cmake_git_revision_result} "unknown")
      endif()
    else()
      set(${vgkit_cmake_git_revision_result} "unknown")
    endif()
  endif()
  unset(vgkit_cmake_git_revision_result)
endmacro()

#
# vgkit_cmake_git_description() defines variable ${PROJECT_NAME_UPPER}_GIT_DESCRIPTION,
#   where ${PROJECT_NAME_UPPER} is the upper-case version of ${PROJECT_NAME}, to
#   the human-readable form produced by `git describe --dirty` command.
#   If the variable is already defined will print a warning message but not change the value of
#   ${PROJECT_NAME_UPPER}_GIT_DESCRIPTION
#
macro(vgkit_cmake_git_description)
  string(TOUPPER "${PROJECT_NAME}_GIT_DESCRIPTION" vgkit_cmake_git_description_result)
  if (DEFINED ${vgkit_cmake_git_description_result})
    message(WARNING "vgkit_cmake_git_description() called for project ${PROJECT_NAME} but the result variable ${vgkit_cmake_git_description_result} already defined; will not change its value")
  else()
    if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
      find_package(Git REQUIRED)
      if (GIT_FOUND)
        execute_process(
                COMMAND ${GIT_EXECUTABLE} describe --dirty
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                RESULT_VARIABLE vgkit_cmake_git_description_errcod
                OUTPUT_VARIABLE ${vgkit_cmake_git_description_result}
                ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT vgkit_cmake_git_description_errcod EQUAL 0)
          set(${vgkit_cmake_git_description_errcod} "unknown")
        endif()
        unset(vgkit_cmake_git_description_errcod)
      else()
        set(${vgkit_cmake_git_description_errcod} "unknown")
      endif()
    else()
      set(${vgkit_cmake_git_description_errcod} "unknown")
    endif()
  endif()
  unset(vgkit_cmake_git_description_errcod)
endmacro()

#
# vgkit_cmake_git_metadata() invokes vgkit_cmake_git_revision and vgkit_cmake_git_description
#
macro(vgkit_cmake_git_metadata)
  vgkit_cmake_git_revision()
  vgkit_cmake_git_description()
endmacro()

cmake_minimum_required(VERSION 3.0.0)

message(STATUS "Resolving GIT Version")

set(_build_version "unknown")


find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE _build_version_only
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} diff --quiet 
    WORKING_DIRECTORY "${local_dir}"
    RESULT_VARIABLE differs
  )
  
# add a trailing plus if the code was change compared to the revision
  if (${differs} EQUAL 0)
#    message(STATUS "GIT_DIFF not differing") 
    string(CONCAT _build_version ${_build_version_only})
  else()
#    message(STATUS "GIT_DIFF differing")
    string(CONCAT _build_version ${_build_version_only}+)
  endif() 
  message(STATUS "GIT hash: ${_build_version}")
else()
  message(STATUS "GIT not found")
endif()

string(TIMESTAMP _time_stamp)

#configure_file(${local_dir}/git_version.h.in ${output_dir}/git_version.h @ONLY)

#configure_file(${local_dir}/gitinfo.cc.in ${output_dir}/gitinfo.cc @ONLY)

set(VERSION "const char* MADNESS_GITREVISION=\"${_build_version}\";")

set(VERSIONFILE "${output_dir}/gitversion.h")
set(VERSIONSOURCE "${source_dir}/gitinfo.cc")

#message(STATUS "looking for file ${VERSIONFILE}")

if(EXISTS ${VERSIONFILE})
#    message(STATUS "reading file ${VERSIONFILE}")
    file(READ ${VERSIONFILE} VERSION_)
    message(STATUS "contents of the gitversion file: ${VERSION_}")
else()
    set(VERSION_ "")
    message(STATUS "cannot find file ${VERSIONFILE}")
endif()


if (NOT "${VERSION}" STREQUAL "${VERSION_}")
    file(WRITE ${VERSIONFILE} "${VERSION}")
    message(STATUS "creating file ${VERSIONFILE}")
    execute_process(
    	COMMAND touch "gitinfo.cc"
    	WORKING_DIRECTORY "${source_dir}"
    )
    
else()
    message(STATUS "file ${FILE} is up to date")
endif()



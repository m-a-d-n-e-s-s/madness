# - CMAKE Config file for the MADNESS package
# The following variables are defined:
#  MADNESS_INCLUDE_DIRS         - The MADNESS include directory
#  MADNESS_LIBRARIES            - The MADNESS libraries and their dependencies
#  MADNESS_<COMPONENT>_FOUND    - System has the specified MADNESS COMPONENT
#  MADNESS_<COMPONENT>_LIBRARY  - The MADNESS COMPONENT library
#  MADNESS_COMPILE_FLAGS        - Compile flags required to build with MADNESS
#  MADNESS_LINKER_FLAGS         - Linker flags required to link with MADNESS
#  MADNESS_VERSION              - MADNESS version number
#  MADNESS_F77_INTEGER_SIZE     - The default F77 integer size used for BLAS calls 

# Compute paths
get_filename_component(MADNESS_CONFIG_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(prefix "/homes/naromero/mpich2-1.5-GCC4.7.3")
set(exec_prefix "${prefix}")
set(MADNESS_DIR "/homes/naromero/mpich2-1.5-GCC4.7.3")
set(MADNESS_INCLUDE_DIRS "${prefix}/include")
set(MADNESS_LIBRARY "${exec_prefix}/lib")

# Set package version
set(MADNESS_VERSION "0.10.1")

# Set compile and link flags, and remove optimization and debug flags
string(REGEX REPLACE "-(O[0-9s]|g[0-9]?)([ ]+|$)" "" MADNESS_COMPILE_FLAGS " -DMPICH_SKIP_MPICXX=1 -DOMPI_SKIP_MPICXX=1  -std=c++11 -O3 -Wall -Wno-strict-aliasing -Wno-deprecated -Wno-unused-local-typedefs  -ffast-math -march=native")
string(REGEX REPLACE "-(O[0-9s]|g[0-9]?)([ ]+|$)" "" MADNESS_LINKER_FLAGS " -std=c++11 -O3 -Wall -Wno-strict-aliasing -Wno-deprecated -Wno-unused-local-typedefs  -ffast-math -march=native ")

# Set MADNESS component variables
set(MADNESS_DEFAULT_COMPONENT_LIST MADchem MADmra MADtinyxml MADmuparser MADlinalg MADtensor MADmisc MADworld)
set(MADNESS_MADchem_DEP_LIST MADmra)
set(MADNESS_MADmra_DEP_LIST MADtinyxml MADmuparser MADlinalg)
set(MADNESS_MADlinalg_DEP_LIST MADtensor)
set(MADNESS_MADtensor_DEP_LIST MADmisc)
set(MADNESS_MADmisc_DEP_LIST MADworld)
set(MADNESS_MADworld_DEP_LIST)

# Check for valid component list
foreach(_comp ${MADNESS_FIND_COMPONENTS})
    if(NOT "${MADNESS_DEFAULT_COMPONENT_LIST}" MATCHES "${_comp}")
        message(FATAL_ERROR "Invalid MADNESS component: ${_comp}")
    endif()
endforeach()


# Set MADNESS libraries variable
foreach(_comp ${MADNESS_DEFAULT_COMPONENT_LIST})

    # Search for MADNESS library
    find_library(MADNESS_${_comp}_LIBRARY ${_comp}
                 PATHS ${MADNESS_LIBRARY}
                 NO_DEFAULT_PATH)

    # Check that the library component was found
    if(MADNESS_${_comp}_LIBRARY)
        set(MADNESS_${_comp}_FOUND TRUE)
    
        # Set MADNESS libraries variable
        if("${MADNESS_FIND_COMPONENTS}" MATCHES "${_comp}" OR NOT MADNESS_FIND_COMPONENTS)
            list(APPEND MADNESS_LIBRARIES ${MADNESS_${_comp}_LIBRARY})
        endif()
    else()
        if(MADNESS_FIND_REQUIRED_${_comp} OR (NOT MADNESS_FIND_COMPONENTS AND MADNESS_FIND_REQUIRED))
            # Fail due to missing required component.
            MESSAGE(FATAL_ERROR "!!ERROR: MADNESS ${_comp} library is not available.")
        endif()
        set(MADNESS_${_comp}_FOUND FALSE)
    endif()

    # Check for dependencies in the component list
    if("${MADNESS_FIND_COMPONENTS}" MATCHES "${_comp}")
        foreach(_comp_dep ${MADNESS_${_comp}_DEP_LIST})
            # Add dependency to the component list if missing
            if(NOT "${MADNESS_FIND_COMPONENTS}" MATCHES "${_comp_dep}")
                list(APPEND MADNESS_FIND_COMPONENTS ${_comp_dep})
            endif()

            # Set required flag for component dependencies
            if(MADNESS_FIND_REQUIRED_${_comp})
                set(MADNESS_FIND_REQUIRED_${_comp_dep} TRUE)
            else()
                set(MADNESS_FIND_REQUIRED_${_comp_dep} FALSE)
            endif()
        endforeach()
    endif()

endforeach()

list(APPEND MADNESS_LIBRARIES " -L/opt/intel/mkl/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -Wl,--end-group -lpthread -lm -ldl  ")

# Set Fortran 77 integer size used by MADNESS
set(MADNESS_F77_INTEGER_SIZE 4)

# Clear local variables
unset(MADNESS_DEFAULT_COMPONENT_LIST)
unset(MADNESS_MADchem_DEP_LIST)
unset(MADNESS_MADmra_DEP_LIST)
unset(MADNESS_MADlinalg_DEP_LIST)
unset(MADNESS_MADtensor_DEP_LIST)
unset(MADNESS_MADmisc_DEP_LIST)
unset(MADNESS_MADworld_DEP_LIST)
unset(prefix)
unset(exec_prefix)

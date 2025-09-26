if (NOT TARGET cereal)
    find_package(cereal QUIET CONFIG)

    if (TARGET cereal)
        message(STATUS "Found Cereal: cereal_CONFIG=${cereal_CONFIG}")
        target_compile_definitions(cereal INTERFACE
                "CEREAL_THREAD_SAFE=1")
    else (TARGET cereal)
        # try looking for header-only cereal
        find_path(cereal_header_paths_tmp
                NAMES
                cereal.hpp
                PATH_SUFFIXES
                include
                cereal
                cereal/include
                )

        get_filename_component(cereal_INCLUDE_DIRS ${cereal_header_paths_tmp} PATH)

        include(FindPackageHandleStandardArgs)
        set(cereal_FIND_QUIETLY 1)
        find_package_handle_standard_args(cereal
                REQUIRED_VARS cereal_INCLUDE_DIRS)

        if (cereal_FOUND)
            message(STATUS "Found Cereal (header-only, at ${cereal_INCLUDE_DIRS})")
            add_library(cereal INTERFACE IMPORTED)
            set_target_properties(cereal PROPERTIES
                    INTERFACE_INCLUDE_DIRECTORIES "${cereal_INCLUDE_DIRS}"
                    INTERFACE_COMPILE_DEFINITIONS "CEREAL_THREAD_SAFE=1")
        endif ()

        # if header-only cereal not found either fetchcontent it over
        if (NOT TARGET cereal)
            cmake_minimum_required(VERSION 3.14.0)  # for FetchContent_MakeAvailable
            include(FetchContent)
            FetchContent_Declare(
                    cereal
                    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
                    GIT_TAG v1.3.2)

            # configure cereal
            set(JUST_INSTALL_CEREAL ON CACHE BOOL "")
            set(THREAD_SAFE ON CACHE BOOL "")
            set(CEREAL_INSTALL ON CACHE BOOL "")

            FetchContent_MakeAvailable(cereal)

            # Export cereal target for build tree so MADNESS can export targets that depend on it
            # Note: cereal v1.3.2 has its own install exports, but we need build tree export for MADNESS
            export(TARGETS cereal FILE "${PROJECT_BINARY_DIR}/cereal-targets.cmake")

            # set cereal_CONFIG to the install location so that we know where to find it
            set(cereal_CONFIG ${CMAKE_INSTALL_PREFIX}/share/cmake/cereal/cereal-config.cmake)
        endif (NOT TARGET cereal)

    endif (TARGET cereal)
endif (NOT TARGET cereal)

if (TARGET cereal)
    set(MADNESS_HAS_CEREAL ON CACHE BOOL "MADNESS has access to Cereal")
else (TARGET cereal)
    message(FATAL_ERROR "MADNESS_ENABLE_CEREAL=ON but could not find or fetch Cereal")
endif (TARGET cereal)
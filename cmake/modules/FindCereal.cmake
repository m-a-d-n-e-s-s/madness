find_package( cereal QUIET )
if( cereal_FOUND )
    target_compile_definitions(cereal INTERFACE
      "CEREAL_THREAD_SAFE=1;JUST_INSTALL_CEREAL=ON;SKIP_PORTABILITY_TEST=ON")
else()
    message( STATUS "Could not find Cereal. Downloading..." )
    include(FetchContent)
    FetchContent_Declare(
      cereal
      GIT_REPOSITORY git@github.com:USCiLab/cereal.git
      GIT_TAG v1.3.0)
    FetchContent_GetProperties(cereal)
    if(NOT cereal_POPULATED)
        FetchContent_Populate( cereal )
        add_library( cereal INTERFACE IMPORTED )
        set_target_properties( cereal PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${cereal_SOURCE_DIR}/include"
          INTERFACE_COMPILE_DEFINITIONS "CEREAL_THREAD_SAFE=1;GAUXC_ENABLE_CEREAL=1;JUST_INSTALL_CEREAL=ON;SKIP_PORTABILITY_TEST=ON")
  endif()
endif()
add_compile_definitions(MADNESS_HAS_CEREAL)


# src

add_subdirectory(madness)
if(ENABLE_APPS)
	add_subdirectory(apps)
endif()
if(ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif()
install(FILES madness.h
    DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}"
    COMPONENT madness)

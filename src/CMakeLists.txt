# src

add_subdirectory(madness)
if(NOT(MADNESS_BUILD_MADWORLD_ONLY OR MADNESS_BUILD_MADWORLD_CLAPACK_ONLY))
	add_subdirectory(apps)
	add_subdirectory(examples)
endif()
install(FILES madness.h
    	DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}"
    	COMPONENT madness)

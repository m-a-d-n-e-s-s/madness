# src

add_subdirectory(madness)
if(NOT (MADNESS_BUILD_MADWORLD_ONLY OR MADNESS_BUILD_LIBRARIES_ONLY) AND LAPACK_FOUND
   AND MADNESS_MRA_HAS_CHEM_DIMENSIONS)
	add_subdirectory(apps)
	add_subdirectory(examples)
endif()
if(MADNESS_HAS_PYTHON)
	if (MADNESS_MRA_HAS_ALL_DIMENSIONS)
		add_subdirectory(pymadness)
	else()
		message(WARNING "pymadness binds all MRA dimensions 1..6; skipped because MADNESS_MRA_DIMENSIONS=${MADNESS_MRA_DIMENSIONS}")
	endif()
endif()
install(FILES madness.h
    	DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}"
    	COMPONENT madness)

if (MADNESS_MRA_HAS_CHEM_DIMENSIONS)
	install(FILES madchem.h
			DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}"
			COMPONENT chem)
endif()

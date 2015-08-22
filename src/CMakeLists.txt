# src

add_subdirectory(madness)
add_subdirectory(apps)
add_subdirectory(examples)

install(FILES madness.h
    DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}"
    COMPONENT madness)

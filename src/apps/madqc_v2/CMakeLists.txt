# src/apps/moldft

set(library_sources
	InputWriter.cpp
)
set(library_headers
	InputWriter.hpp
)

set(MY_EXECUTABLES
	madqc2
)


add_mad_library(qcapp library_sources library_headers "chem" "molresponse2/")
set(targetname qcapp)

foreach (EXECUTABLE ${MY_EXECUTABLES})
    add_mad_executable(${EXECUTABLE} "${EXECUTABLE}.cpp" "MADresponse2;MADchem;MADqcapp")
    target_compile_features(${EXECUTABLE} PRIVATE cxx_std_17)
endforeach ()


install(TARGETS madqc2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
# For testing script




# output directory
file(MAKE_DIRECTORY  ${PacBioConsensus_SwigDir}/lib)

configure_file(
    ${PacBioConsensus_RootDir}/swig/ConsensusCore2.i.in
    ${PacBioConsensus_RootDir}/swig/ConsensusCore2.i
)

include_directories(
    ${PacBioConsensus_IncludeDir}
)

include(FindPythonPaths)
find_python_inc(PYTHON_INCLUDE_DIRS)
find_numpy_inc(NUMPY_INCLUDE_DIRS)

# Generate paths for test data
#configure_file(
#    ${PacBioConsensus_TestsDir}/TestData.h.in
#    ${PacBioConsensus_TestsDir}/TestData.h
#)

include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    ${NUMPY_INCLUDE_DIRS}
)

# fPIC is required by a the majority of Linux systems,
# and the original ConsensusCore used it by default.
# We could use python -c "print(sysconfig.get_config_var('CCSHARED'))",
# but it's honestly unnecessary and requires providing PYTHON_EXECUTABLE
target_compile_options(pbconsensus PUBLIC -fPIC)
add_compile_options(-fPIC)

# ConsensusCore2_wrap.cxx
set(SWIG_COMMAND swig CACHE PATH "swig path")

add_custom_command(OUTPUT ConsensusCore2_wrap.cxx
    COMMAND ${SWIG_COMMAND} -Wextra -c++ -python -builtin -module ConsensusCore2
        -I${PacBioConsensus_IncludeDir}
        -o ${CMAKE_BINARY_DIR}/swig/ConsensusCore2_wrap.cxx -outdir ${CMAKE_BINARY_DIR}/swig
        ${PacBioConsensus_SwigDir}/ConsensusCore2.i
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/swig/ConsensusCore2.py
        ${PacBioConsensus_SwigDir}/lib/__init__.py
)

# _ConsensusCore2.so
add_library(_ConsensusCore2 MODULE
    ConsensusCore2_wrap.cxx
)

target_link_libraries(_ConsensusCore2
    ${PBCONSENSUS_LIBRARIES}
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(_ConsensusCore2 "-Wl,-undefined,dynamic_lookup")
endif()

set_target_properties(_ConsensusCore2 PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PacBioConsensus_SwigDir}/lib
    PREFIX ""
)

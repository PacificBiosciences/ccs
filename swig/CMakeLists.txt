
# output directory
file(MAKE_DIRECTORY  ${PacBioConsensus_SwigDir}/lib)

include_directories(
    ${PacBioConsensus_IncludeDir}
)

if (NOT PYTHON_INCLUDE_DIRS OR
    NOT PYTHON_LIBRARIES)
    find_package(PythonLibs REQUIRED)
endif()

if (NOT NUMPY_INCLUDE_DIRS)
    message(FATAL_ERROR "Need to provide NUMPY_INCLUDE_DIRS")
endif()

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioConsensus_CXX_FLAGS}")

# pbconsensus_wrap.cxx
add_custom_command(OUTPUT pbconsensus_wrap.cxx
    COMMAND swig -Wall -c++ -python -builtin -module pbconsensus
                 -I${PacBioConsensus_IncludeDir}
                 -o ${CMAKE_BINARY_DIR}/swig/pbconsensus_wrap.cxx -outdir ${CMAKE_BINARY_DIR}/swig
                 ${PacBioConsensus_SwigDir}/pbconsensus.i
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/swig/pbconsensus.py
                                     ${PacBioConsensus_SwigDir}/lib
)

# _pbconsensus.so
add_library(_pbconsensus MODULE
    pbconsensus_wrap.cxx
)

target_link_libraries(_pbconsensus
    ${PBCONSENSUS_STATIC}
    ${PYTHON_LIBRARIES}
)

set_target_properties(_pbconsensus PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PacBioConsensus_SwigDir}/lib
    PREFIX ""
)


# output directory
file(MAKE_DIRECTORY  ${PacBioConsensus_SwigDir}/lib)

include_directories(
    ${PacBioConsensus_IncludeDir}
)

if (NOT PYTHON_INCLUDE_DIRS OR
    NOT PYTHON_LIBRARIES)
    if (NOT PYTHON_EXECUTABLE)
        find_package(PythonInterp REQUIRED)
    endif()
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_inc(), end='')"
        RESULT_VARIABLE PYTHON_INCLUDE_SUCCESS
        OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import distutils.sysconfig; from os.path import dirname as d; print(d(distutils.sysconfig.get_python_lib(standard_lib=True)), end='')"
        RESULT_VARIABLE PYTHON_LIBPATH_SUCCESS
        OUTPUT_VARIABLE PYTHON_LIBPATH)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_version(), end='')"
        RESULT_VARIABLE PYTHON_VERSION_SUCCESS
        OUTPUT_VARIABLE PYTHON_VERSION)
    # find the library (first static, then shared, then module)
    set(__PYTHON_LIBRARY_PREFIX__ "${PYTHON_LIBPATH}/libpython${PYTHON_VERSION}")
    file(GLOB PYTHON_LIBRARIES "${__PYTHON_LIBRARY_PREFIX__}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (NOT PYTHON_LIBRARIES)
        file(GLOB PYTHON_LIBRARIES "${__PYTHON_LIBRARY_PREFIX__}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
    if (NOT PYTHON_LIBRARIES)
        file(GLOB PYTHON_LIBRARIES "${__PYTHON_LIBRARY_PREFIX__}${CMAKE_SHARED_MODULE_SUFFIX}")
    endif()
    unset(__PTYHON_LIBRARY_PREFIX__)
    # check for success
    if (NOT PYTHON_INCLUDE_SUCCESS MATCHES 0 OR
        NOT PYTHON_LIBPATH_SUCCESS MATCHES 0 OR
        NOT PYTHON_VERSION_SUCCESS MATCHES 0 OR
        NOT PYTHON_INCLUDE_DIRS OR
        NOT PYTHON_LIBRARIES)
        message(FATAL_ERROR "PYTHON_INCLUDE_DIRS and PYTHON_LIBRARIES need to be set manually")
    endif()
    message(STATUS "Found PythonLibs: ${PYTHON_LIBRARIES}")
endif()

if (NOT NUMPY_INCLUDE_DIRS)
    if (NOT PYTHON_EXECUTABLE)
        find_package(PythonInterp REQUIRED)
    endif()
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import numpy; print(numpy.get_include(), end='')"
        RESULT_VARIABLE NUMPY_INCLUDE_SUCCESS
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS)
    if (NOT NUMPY_INCLUDE_SUCCESS MATCHES 0 OR
        NOT NUMPY_INCLUDE_DIRS)
        message(FATAL_ERROR "NUMPY_INCLUDE_DIRS needs to be set manually")
    endif()
    message(STATUS "Found Numpy: ${NUMPY_INCLUDE_DIRS}")
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
                                     ${PacBioConsensus_SwigDir}/lib/__init__.py
)

# _pbconsensus.so
add_library(_pbconsensus MODULE EXCLUDE_FROM_ALL
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

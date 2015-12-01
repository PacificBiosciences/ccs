
if (__find_python_inc_lib)
    return()
endif()
set(__find_python_inc_lib YES)

function(find_python_inc_lib _PYTHON_INC _PYTHON_LIB)
    # find the executable
    if (NOT PYTHON_EXECUTABLE)
        find_package(PythonInterp REQUIRED)
    endif()
    # find the include directory
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_inc(), end='')"
        RESULT_VARIABLE PYTHON_INCLUDE_SUCCESS
        OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS)
    # find potential library paths (curse you Debian!)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function
from distutils import sysconfig as sc
from os.path import dirname as d
print(sc.get_config_var('LIBDIR') + ';' + d(sc.get_config_var('LIBPC')), end='')"
        RESULT_VARIABLE PYTHON_LIBPATH_SUCCESS
        OUTPUT_VARIABLE PYTHON_LIBPATHS)
    # find the python version
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_version(), end='')"
        RESULT_VARIABLE PYTHON_VERSION_SUCCESS
        OUTPUT_VARIABLE PYTHON_VERSION)
    # find the library (first shared, then module)
    set(__PYTHON_LIBRARY_SUFFIXES__ ${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_SHARED_MODULE_SUFFIX})
    foreach(PYTHON_LIBRARY_SUFFIX IN LISTS __PYTHON_LIBRARY_SUFFIXES__)
        foreach(PYTHON_LIBPATH IN LISTS PYTHON_LIBPATHS)
            file(GLOB PYTHON_LIBRARIES "${PYTHON_LIBPATH}/libpython${PYTHON_VERSION}${PYTHON_LIBRARY_SUFFIX}")
            if (PYTHON_LIBRARIES)
                break()
            endif()
        endforeach()
        if (PYTHON_LIBRARIES)
            break()
        endif()
    endforeach()
    # check for success
    if (NOT PYTHON_INCLUDE_SUCCESS EQUAL 0 OR
        NOT PYTHON_LIBPATH_SUCCESS EQUAL 0 OR
        NOT PYTHON_VERSION_SUCCESS EQUAL 0 OR
        NOT PYTHON_INCLUDE_DIRS OR
        NOT PYTHON_LIBRARIES)
        message(FATAL_ERROR "${_PYTHON_INC} and ${_PYTHON_LIB} need to be set manually")
    endif()
    # set the output variables
    set(${_PYTHON_INC} "${PYTHON_INCLUDE_DIRS}" PARENT_SCOPE)
    set(${_PYTHON_LIB} "${PYTHON_LIBRARIES}" PARENT_SCOPE)
endfunction()

function(find_numpy_inc _NUMPY_INC)
    if (NOT PYTHON_EXECUTABLE)
        find_package(PythonInterp REQUIRED)
    endif()
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import numpy; print(numpy.get_include(), end='')"
        RESULT_VARIABLE NUMPY_INCLUDE_SUCCESS
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS)
    if (NOT NUMPY_INCLUDE_SUCCESS EQUAL 0 OR
        NOT NUMPY_INCLUDE_DIRS)
        message(FATAL_ERROR "NUMPY_INCLUDE_DIRS needs to be set manually")
    endif()
    set(${_NUMPY_INC} "${NUMPY_INCLUDE_DIRS}" PARENT_SCOPE)
endfunction()


if (__find_python_inc)
    return()
endif()
set(__find_python_inc YES)

function(find_python_inc _PYTHON_INC)
    # find the executable
    if (NOT PYTHON_EXECUTABLE)
        find_package(PythonInterp REQUIRED)
    endif()
    # find the include directory
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
        "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_inc(), end='')"
        RESULT_VARIABLE PYTHON_INCLUDE_SUCCESS
        OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS)
    # check for success
    if (NOT PYTHON_INCLUDE_SUCCESS EQUAL 0 OR
        NOT PYTHON_INCLUDE_DIRS)
        message(FATAL_ERROR "${_PYTHON_INC} needs to be set manually")
    endif()
    # set the output variables
    set(${_PYTHON_INC} "${PYTHON_INCLUDE_DIRS}" PARENT_SCOPE)
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

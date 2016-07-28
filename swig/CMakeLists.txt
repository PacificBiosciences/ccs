
# output directory
file(MAKE_DIRECTORY  ${CMAKE_BINARY_DIR}/swig/lib)

include(uny-findpythonpath)
find_python_inc(PYTHON_INCLUDE_DIRS)
find_numpy_inc(NUMPY_INCLUDE_DIRS)

include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
    ${NUMPY_INCLUDE_DIRS}
)

# fPIC is required by a the majority of Linux systems,
# and the original ConsensusCore used it by default.
# We could use python -c "print(sysconfig.get_config_var('CCSHARED'))",
# but it's honestly unnecessary and requires providing PYTHON_EXECUTABLE
target_compile_options(cc2 PUBLIC -fPIC)
add_compile_options(-fPIC)

# ConsensusCore2_wrap.cxx
set(SWIG_COMMAND swig CACHE PATH "swig PATH")

add_custom_command(OUTPUT ConsensusCore2_wrap.cxx
    COMMAND ${SWIG_COMMAND} -Wextra -c++ -python -builtin -module ConsensusCore2
        -I${UNY_IncludeDir} -I${CMAKE_BINARY_DIR}/generated
        -o ${CMAKE_BINARY_DIR}/swig/ConsensusCore2_wrap.cxx -outdir ${CMAKE_BINARY_DIR}/swig
        ${UNY_SwigDir}/ConsensusCore2.i
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/swig/ConsensusCore2.py
        ${CMAKE_BINARY_DIR}/swig/lib/__init__.py
)

# _ConsensusCore2.so
add_library(_ConsensusCore2 MODULE
    ConsensusCore2_wrap.cxx
)

target_include_directories(_ConsensusCore2
    PUBLIC ${UNY_IncludeDir}
    PUBLIC ${CMAKE_BINARY_DIR}/generated
)

set(CC2_LIBRARIES
    ${PRE_LINK}
    cc2
    ${POST_LINK}
    CACHE INTERNAL "" FORCE)

target_link_libraries(_ConsensusCore2
    ${CC2_LIBRARIES}
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(_ConsensusCore2 "-Wl,-undefined,dynamic_lookup")
endif()

set_target_properties(_ConsensusCore2 PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/swig/lib
    PREFIX ""
)

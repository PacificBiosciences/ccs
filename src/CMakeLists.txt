# get all header files for IDE support
file(GLOB_RECURSE UNY_HEADER "${UNY_IncludeDir}/*.h")
file(GLOB_RECURSE UNY_HIDDEN_HEADER "*.h")
# get all cpp files
file(GLOB_RECURSE UNY_CPP "*.cpp")
# exclude executable cpp fles
file(GLOB_RECURSE UNY_MAIN_CPP "main/*.cpp")
list(REMOVE_ITEM UNY_CPP ${UNY_MAIN_CPP})

# add main library including everything
add_library(unanimity STATIC
    ${UNY_CPP}
    ${UNY_HEADER}
    ${UNY_HIDDEN_HEADER}
    ${CPPOPTPARSE_CPP}
)

# includes
target_include_directories(unanimity PUBLIC
    ${UNY_IncludeDir}
    ${CMAKE_BINARY_DIR}/generated
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${UNY_ThirdPartyDir}/seqan/include
    ${CPPOPTPARSE_H}
    ${PacBioBAM_INCLUDE_DIRS}
)

if(APPLE)
    set(PRE_LINK -force_load)
elseif(UNIX)
    SET(PRE_LINK -Wl,-whole-archive)
    SET(POST_LINK -Wl,-no-whole-archive)
endif()

# because we can be a dependency, set this in the parent scope if so
if(${ROOT_PROJECT_NAME} STREQUAL "UNANIMITY")
    set(UNANIMITY_LIBRARIES ${PRE_LINK} unanimity ${POST_LINK} CACHE INTERNAL "" FORCE)
else()
    set(UNANIMITY_LIBRARIES ${PRE_LINK} unanimity ${POST_LINK} PARENT_SCOPE CACHE INTERNAL "" FORCE)
endif()

# add executables
if (UNY_build_bin)
    add_executable(ccs ${UNY_SourceDir}/main/ccs.cpp)
    
    set_target_properties(ccs PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    target_link_libraries(ccs 
        ${UNANIMITY_LIBRARIES}
        ${PacBioBAM_LIBRARIES}
        ${HTSLIB_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        ${CMAKE_DL_LIBS}
        ${ZLIB_LIBRARIES}
    )

    if (TARGET htslib)
        add_dependencies(ccs htslib)
    endif()

    install(TARGETS ccs RUNTIME DESTINATION bin)
endif()
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
set(UNY_INCLUDE_DIRS
    ${UNY_IncludeDir}
    ${CMAKE_BINARY_DIR}/generated
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${SEQAN_INCLUDE_DIRS}
    ${CPPOPTPARSE_IncludeDir}
    ${PacBioBAM_INCLUDE_DIRS}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

target_include_directories(unanimity PUBLIC
    ${UNY_INCLUDE_DIRS}
)

if(APPLE)
    set(PRE_LINK -force_load)
elseif(UNIX)
    SET(PRE_LINK -Wl,-whole-archive)
    SET(POST_LINK -Wl,-no-whole-archive)
endif()

set(UNY_LIBRARIES 
    ${PRE_LINK} 
    unanimity
    ${POST_LINK}
    ${ZLIB_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    ${PacBioBAM_LIBRARIES}
    CACHE INTERNAL "" FORCE)

if (TARGET htslib)
    add_dependencies(unanimity htslib)
endif()

if (TARGET pbbam)
    add_dependencies(unanimity pbbam)
endif()

# add executables
if (UNY_build_bin)
    add_executable(ccs ${UNY_SourceDir}/main/ccs.cpp)
    
    set_target_properties(ccs PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    target_link_libraries(ccs 
        ${UNY_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        ${CMAKE_DL_LIBS}
    )

    if (TARGET htslib)
        add_dependencies(ccs htslib)
    endif()

    if (TARGET pbbam)
        add_dependencies(ccs pbbam)
    endif()

    install(TARGETS ccs RUNTIME DESTINATION bin)
endif()
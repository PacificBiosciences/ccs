# get all cpp files
file(GLOB_RECURSE UNY_CPP "*.cpp")
# exclude executable cpp fles
file(GLOB_RECURSE UNY_MAIN_CPP "main/*.cpp")
list(REMOVE_ITEM UNY_CPP ${UNY_MAIN_CPP})

# add main library including everything
add_library(unanimity STATIC
    ${UNY_CPP}
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
)

# link dependency libs
target_link_libraries(unanimity
    pbbam
    pbcopper
    ${HTSLIB_LIBRARIES}
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
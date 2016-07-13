##############################################
# CMake build script for the UNANIMITY library
##############################################

cmake_policy(SET CMP0048 NEW)
project(UNANIMITY VERSION 0.0.1 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.2)

set(PacBioBAM_build_docs    OFF CACHE INTERNAL "" FORCE)
set(PacBioBAM_build_tests   OFF CACHE INTERNAL "" FORCE)
set(PacBioBAM_build_tools   OFF CACHE INTERNAL "" FORCE)

# Build type
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile RelWithDebInfo" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Build-time options
option(UNY_build_bin    "Build binaries." ON)
option(UNY_build_tests  "Build UNANIMITY's unit tests." ON)
option(UNY_use_ccache   "Build UNANIMITY using ccache, if available." ON)

# Main project paths
set(UNY_RootDir       ${UNANIMITY_SOURCE_DIR})
set(UNY_IncludeDir    ${UNY_RootDir}/include)
set(UNY_SourceDir     ${UNY_RootDir}/src)
set(UNY_TestsDir      ${UNY_RootDir}/test)
set(UNY_ThirdPartyDir ${UNY_RootDir}/third-party)

# Project configuration
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})
include(uny-ccache)
include(uny-compilerflags)
include(uny-dependencies)
include(uny-gitsha1)

# Config generation
find_git_sha1(UNY_GIT_SHA1)
file (STRINGS "CHANGELOG.md" UNY_CHANGELOG)
configure_file(
    ${UNY_IncludeDir}/pacbio/unanimity/Version.h.in
    ${CMAKE_BINARY_DIR}/generated/pacbio/unanimity/Version.h
)

# shared CXX flags for src & tests
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UNY_FLAGS}")

# Build library
# add_subdirectory(src/consensus)

# Build tests
# if(UNY_build_tests)
#     add_subdirectory(test/unit)
# endif()

# Build binaries
# if (UNY_build_bin)
#     add_subdirectory(src/pbccs)

#     install(TARGETS ${ALL_EXE_TARGETS} RUNTIME DESTINATION bin)
# endif()
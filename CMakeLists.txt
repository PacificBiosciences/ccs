########################################################################
# CMake build script for PacBioConsensus library.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set the version in project()
project(PacBioConsensus VERSION 0.12.0 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.0)

# set magic variable
set(ROOT_PROJECT_NAME ${PROJECT_NAME} CACHE STRING "root project name")

# required packages
if (NOT Boost_INCLUDE_DIRS)
    find_package(Boost 1.58.0 REQUIRED)
endif()

# main project paths
set(PacBioConsensus_RootDir           ${PacBioConsensus_SOURCE_DIR})
set(PacBioConsensus_IncludeDir        ${PacBioConsensus_RootDir}/include)
set(PacBioConsensus_LibDir            ${PacBioConsensus_RootDir}/lib)
set(PacBioConsensus_SourceDir         ${PacBioConsensus_RootDir}/src)
set(PacBioConsensus_SwigDir           ${PacBioConsensus_RootDir}/swig)
set(PacBioConsensus_TestsDir          ${PacBioConsensus_RootDir}/tests)

file(MAKE_DIRECTORY                   ${PacBioConsensus_LibDir})

# find the git sha1
list(APPEND CMAKE_MODULE_PATH "${PacBioConsensus_RootDir}/cmake")
include(FindGitSha1)
find_git_sha1(CC2_GIT_SHA1)

configure_file(
    ${PacBioConsensus_IncludeDir}/pacbio/consensus/Version.h.in
    ${PacBioConsensus_IncludeDir}/pacbio/consensus/Version.h
)

file(GLOB_RECURSE PacBioConsensus_CPP "${PacBioConsensus_SourceDir}/*.cpp")
file(GLOB_RECURSE PacBioConsensus_HPP "${PacBioConsensus_IncludeDir}/*.h")

# shared CXX flags for src & tests
include(CheckCXXCompilerFlag)
set(PacBioConsensus_CXX_FLAGS "-std=c++11 -Wall -pedantic -msse3")

# silence gmock warning
check_cxx_compiler_flag("-Wno-unused-variable" HAS_NO_UNUSED_VARIABLE)
if(HAS_NO_UNUSED_VARIABLE)
    set(PacBioConsensus_CXX_FLAGS "${PacBioConsensus_CXX_FLAGS} -Wno-unused-variable")
endif(HAS_NO_UNUSED_VARIABLE)

# includes
include_directories(. ${PacBioConsensus_IncludeDir})
include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
)

# make release with debug the default, helps tests pass by being fast enough
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build,
      options are: Debug Release Profile RelWithDebInfo" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioConsensus_CXX_FLAGS}")

# libpbconsensus.a
add_library(pbconsensus STATIC
    ${PacBioConsensus_CPP}
    ${PacBioConsensus_HPP}
)
set_target_properties(pbconsensus PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PacBioConsensus_LibDir}
)

# support static registration
if(APPLE)
    set(PRE_LINK -force_load)
elseif(UNIX)
    SET(PRE_LINK -Wl,-whole-archive)
    SET(POST_LINK -Wl,-no-whole-archive)
endif()

# because we can be a dependency, set this in the parent scope if so
if(${ROOT_PROJECT_NAME} STREQUAL "PacBioConsensus")
    set(PBCONSENSUS_LIBRARIES ${PRE_LINK} pbconsensus ${POST_LINK})
else()
    set(PBCONSENSUS_LIBRARIES ${PRE_LINK} pbconsensus ${POST_LINK} PARENT_SCOPE)
endif()

# swig
if (PYTHON_SWIG)
    add_subdirectory(${PacBioConsensus_SwigDir})
endif()

# testing
add_subdirectory(${PacBioConsensus_TestsDir})

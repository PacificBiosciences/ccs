########################################################################
# CMake build script for the pbchimera library.
########################################################################

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 OLD)
endif()

project(PacBioChimera CXX C)
cmake_minimum_required(VERSION 2.8)

#set magic variable
set(ROOT_PROJECT_NAME ${PROJECT_NAME} CACHE STRING "root project name")

# project version
set(PacBioChimera_VERSION_MAJOR 0)
set(PacBioChimera_VERSION_MINOR 0)
set(PacBioChimera_VERSION_PATCH 1)
set(PacBioChimera_VERSION
  "${PacBioChimera_VERSION_MAJOR}.${PacBioChimera_VERSION_MINOR}.${PacBioChimera_VERSION_PATCH}"
)

# build-time options
#option(PacBioChimera_build_docs  "Build pbchimera's API documentation." ON)
#option(PacBioChimera_build_tests "Build pbchimera's unit tests." ON)

# main project paths
set(PacBioChimera_RootDir       ${PacBioChimera_SOURCE_DIR})
set(PacBioChimera_BinDir        ${PacBioChimera_RootDir}/bin)
set(PacBioChimera_DocsDir       ${PacBioChimera_RootDir}/docs)
set(PacBioChimera_LibDir        ${PacBioChimera_RootDir}/lib)
set(PacBioChimera_IncludeDir    ${PacBioChimera_RootDir}/include)
set(PacBioChimera_SourceDir     ${PacBioChimera_RootDir}/src)
set(PacBioChimera_TestsDir      ${PacBioChimera_RootDir}/test)

# External include directories
set(SEQAN_IncludeDir            ${PacBioChimera_RootDir}/../seqan/include CACHE PATH "Seqan directory")

# set this compiler definition to prevent unhygenic seqan behavior
add_definitions(-DSEQAN_NO_GLOBAL_EXCEPTION_HANDLER)

# shared CXX flags for src & tests
include(CheckCXXCompilerFlag)
set(PacBioChimera_CXX_FLAGS "-std=c++11 -Wall -Wextra -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioChimera_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "-lz")

include_directories(
    ${SEQAN_IncludeDir}
    ${PacBioChimera_IncludeDir}
)

# define symbols for projects that use pbchimera
set(PacBioChimera_INCLUDE_DIRS
    ${PacBioChimera_IncludeDir}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

find_package(Threads REQUIRED)

# Add the sub-directory with
add_subdirectory(src)

# gmock/gtest
if (GMOCK_RootDir)
    add_subdirectory(test/unit)
endif()

# find the git sha1
list(APPEND CMAKE_MODULE_PATH "${PacBioChimera_RootDir}/cmake")
include(FindGitSha1)
find_git_sha1(CHIMERA_GIT_SHA1)
configure_file("${PacBioChimera_IncludeDir}/pacbio/chimera/Version.h.in" "${PacBioChimera_IncludeDir}/pacbio/chimera/Version.h" @ONLY)

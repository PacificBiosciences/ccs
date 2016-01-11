########################################################################
# CMake build script for the pbchimera library.
########################################################################

project(PacBioChimera CXX C)
cmake_minimum_required(VERSION 2.8)

# project version
set(PacBioChimera_MAJOR_VERSION 0)
set(PacBioChimera_MINOR_VERSION 0)
set(PacBioChimera_PATCH_VERSION 1)
set(PacBioChimera_VERSION
  "${PacBioChimera_MAJOR_VERSION}.${PacBioChimera_MINOR_VERSION}.${PacBioChimera_PATCH_VERSION}"
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
set(PacBioChimera_TestsDir      ${PacBioChimera_RootDir}/tests)

# External include directories
set(SEQAN_IncludeDir            ${PacBioChimera_RootDir}/../seqan/include)

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

add_subdirectory(src)

########################################################################
# CMake build script for the pbchimera library.
########################################################################

project(pbchimera CXX C)
cmake_minimum_required(VERSION 2.8)

# project version
set(pbchimera_MAJOR_VERSION 0)
set(pbchimera_MINOR_VERSION 0)
set(pbchimera_PATCH_VERSION 1)
set(pbchimera_VERSION
  "${pbchimera_MAJOR_VERSION}.${pbchimera_MINOR_VERSION}.${pbchimera_PATCH_VERSION}"
)

# build-time options
#option(pbchimera_build_docs  "Build pbchimera's API documentation." ON)
#option(pbchimera_build_tests "Build pbchimera's unit tests." ON)

# main project paths
set(pbchimera_RootDir       ${pbchimera_SOURCE_DIR})
set(pbchimera_BinDir        ${pbchimera_RootDir}/bin)
set(pbchimera_DocsDir       ${pbchimera_RootDir}/docs)
set(pbchimera_LibDir        ${pbchimera_RootDir}/lib)
set(pbchimera_IncludeDir    ${pbchimera_RootDir}/include)
set(pbchimera_SourceDir     ${pbchimera_RootDir}/src)
set(pbchimera_TestsDir      ${pbchimera_RootDir}/tests)
#file(MAKE_DIRECTORY ${pbchimera_LibDir})

find_package( ZLIB REQUIRED )

# Include directories
set(SEQAN_IncludeDir        ${pbchimera_RootDir}/../seqan/include)

# shared CXX flags for src & tests
include(CheckCXXCompilerFlag)
set(pbchimera_CXX_FLAGS "-std=c++11 -Wall -Wextra -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${pbchimera_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "-lz")

# set this compiler definition to prevent unhygenic seqan behavior
add_definitions(-DSEQAN_NO_GLOBAL_EXCEPTION_HANDLER)

include_directories(
    ${ZLIB_INCLUDE_DIRS}
    ${SEQAN_IncludeDir}
    ${pbchimera_IncludeDir}
)

# define symbols for projects that use pbchimera
set(pbchimera_INCLUDE_DIRS
    ${pbchimera_IncludeDir}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

# set up library include dirs
include_directories(
    ${pbchimera_IncludeDir}
)

add_subdirectory(src)

#find_package(Boost 1.36 COMPONENTS unit_test_framework REQUIRED)
#enable_testing()
#add_subdirectory(test/unit)
#add_test(NAME Adapter COMMAND AdapterTest)

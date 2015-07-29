########################################################################
# CMake build script for ConsensusCore library.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set the version in project()
project(ConsensusCore VERSION 1.0.1 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.0)

# required packages
if (NOT Boost_INCLUDE_DIRS)
    find_package(Boost 1.58.0 REQUIRED)
endif()

find_package(Threads REQUIRED)

if (NOT ZLIB_INCLUDE_DIRS OR
    NOT ZLIB_LIBRARIES)
    find_package(ZLIB REQUIRED)
endif()

# project version
set(ConsensusCore_VERSION
  "${ConsensusCore_VERSION_MAJOR}.${ConsensusCore_VERSION_MINOR}.${ConsensusCore_VERSION_PATCH}"
)

# main project paths
set(ConsensusCore_RootDir           ${ConsensusCore_SOURCE_DIR})
set(ConsensusCore_IncludeDir        ${ConsensusCore_RootDir}/include)
set(ConsensusCore_LibDir            ${ConsensusCore_RootDir}/lib)
set(ConsensusCore_SourceDir         ${ConsensusCore_RootDir}/src/C++)
set(ConsensusCore_TestsDir          ${ConsensusCore_RootDir}/src/Tests)
set(ConsensusCore_ThirdPartyDir     ${ConsensusCore_RootDir}/external)

file(MAKE_DIRECTORY                 ${ConsensusCore_LibDir})

file(GLOB_RECURSE ConsensusCore_CPP "${ConsensusCore_RootDir}/src/C++/*.cpp")
file(GLOB_RECURSE ConsensusCore_HPP "${ConsensusCore_RootDir}/include/*.hpp")

# third-party libraries
set(CPPOPTPARSE_RootDir         ${ConsensusCore_ThirdPartyDir}/cpp-optparse)
set(CPPOPTPARSE_IncludeDir      ${CPPOPTPARSE_RootDir})
set(SEQAN_IncludeDir            ${ConsensusCore_ThirdPartyDir}/seqan-library-2.0.0/include)

# shared CXX flags for src & tests
include(CheckCXXCompilerFlag)
set(ConsensusCore_CXX_FLAGS "-std=c++11 -Wall -pedantic -msse3")

# silence gmock warning
check_cxx_compiler_flag("-Wno-unused-variable" HAS_NO_UNUSED_VARIABLE)
if(HAS_NO_UNUSED_VARIABLE)
    set(ConsensusCore_CXX_FLAGS "${ConsensusCore_CXX_FLAGS} -Wno-unused-variable")
endif(HAS_NO_UNUSED_VARIABLE)

# includes
include_directories(. ${ConsensusCore_IncludeDir})
include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${ConsensusCore_ThirdPartyDir}
)

# flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ConsensusCore_CXX_FLAGS}")

# gperftools support
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_pie")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

# libpbccs.a
add_library(ConsensusCore STATIC ${ConsensusCore_CPP} ${ConsensusCore_HPP})
set_target_properties(ConsensusCore PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${ConsensusCore_LibDir}
)

# testing
#add_subdirectory(${ConsensusCore_TestsDir})

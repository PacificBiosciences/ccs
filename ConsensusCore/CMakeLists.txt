########################################################################
# CMake build script for ConsensusCore library.
########################################################################

project(ConsensusCore CXX C)
cmake_minimum_required(VERSION 2.8)

# required packages
find_package(Boost 1.58.0 REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# project version
set(ConsensusCore_MAJOR_VERSION 1)
set(ConsensusCore_MINOR_VERSION 0)
set(ConsensusCore_PATCH_VERSION 1)
set(ConsensusCore_VERSION
  "${ConsensusCore_MAJOR_VERSION}.${ConsensusCore_MINOR_VERSION}.${ConsensusCore_PATCH_VERSION}"
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

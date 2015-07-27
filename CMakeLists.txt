########################################################################
# CMake build script for PacBioCCS executable.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set the version in project()
project(PacBioCCS VERSION 0.0.1 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.0)

# required packages
find_package(Boost 1.58.0 REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# project version
set(PacBioCCS_VERSION
    "${PacBioCCS_VERSION_MAJOR}.${PacBioCCS_VERSION_MINOR}.${PacBioCCS_VERSION_PATCH}"
)

# main project paths
set(PacBioCCS_RootDir           ${PacBioCCS_SOURCE_DIR})
set(PacBioCCS_IncludeDir        ${PacBioCCS_RootDir}/include)
set(PacBioCCS_DocsDir           ${PacBioCCS_RootDir}/docs)
set(PacBioCCS_SourceDir         ${PacBioCCS_RootDir}/src)
set(PacBioCCS_TestsDir          ${PacBioCCS_RootDir}/tests)
set(PacBioCCS_ThirdPartyDir     ${PacBioCCS_RootDir}/third-party)

if(NOT PacBioCCS_OutputDir)
    set(PacBioCCS_OutputDir ${PacBioCCS_RootDir})
endif()

set(PacBioCCS_BinDir            ${PacBioCCS_OutputDir}/bin)
set(PacBioCCS_LibDir            ${PacBioCCS_OutputDir}/lib)

file(MAKE_DIRECTORY             ${PacBioCCS_BinDir})
file(MAKE_DIRECTORY             ${PacBioCCS_LibDir})

# PacBio libraries
set(Bifx_RootDir                ${PacBioCCS_RootDir}/../../../software/smrtanalysis/bioinformatics CACHE PATH "smrtanalysis/bioinformatics root path")

# ConsensusCore
#set(ConsensusCore_RootDir       ${Bifx_RootDir}/lib/cpp/ConsensusCore                              CACHE PATH "ConsensusCore root path")
set(ConsensusCore_RootDir       ${PacBioCCS_RootDir}/ConsensusCore)
set(ConsensusCore_IncludeDir    ${ConsensusCore_RootDir}/include ${ConsensusCore_RootDir}/external)

# pbbam
set(PBBam_RootDir               ${Bifx_RootDir}/staging/PostPrimary/pbbam                          CACHE PATH "pbbam source path")
set(PBBam_IncludeDir            ${PBBam_RootDir}/include)
#set(PBBam_LibDir                ${PBBam_RootDir}/lib)

# htslib
set(HTS_RootDir                 ${Bifx_RootDir}/staging/PostPrimary/htslib)
set(HTS_IncludeDir              ${HTS_RootDir})
#set(HTS_LibDir                  ${HTS_RootDir})

# prebuilt libraries
set(PreBuilt_OutDir             ${Bifx_RootDir}/../prebuilt.out                                    CACHE PATH "smrtanalysis/prebuilt.out root path")

# third-party libraries
set(CPPOPTPARSE_RootDir         ${PacBioCCS_ThirdPartyDir}/cpp-optparse)
set(CPPOPTPARSE_IncludeDir      ${CPPOPTPARSE_RootDir})
set(SEQAN_IncludeDir            ${PacBioCCS_ThirdPartyDir}/seqan-library-2.0.0/include)

# laa sources
file(GLOB         PacBioCCS_CPP "${PacBioCCS_SourceDir}/*.cpp")
file(GLOB_RECURSE PacBioCCS_H   "${PacBioCCS_IncludeDir}/*.h")

# cpp-optparse sources
set(CPPOPTPARSE_CPP             ${CPPOPTPARSE_RootDir}/OptionParser.cpp)

# find (existing) libraries needed by executable and tests
#find_library(CONSENSUSCORE_LIBRARIES
#    NAMES ConsensusCore
#    PATHS ${ConsensusCore_LibDir}
#    NO_DEFAULT_PATH)

#find_library(PBBAM_LIBRARY
#    NAMES pbbam
#    PATHS ${PBBam_LibDir} ${HTS_LibDir}
#    NO_DEFAULT_PATH)
#
#find_library(HTS_LIBRARY
#    NAMES hts
#    PATHS ${HTS_LibDir}
#    NO_DEFAULT_PATH)

#set(PBBAM_LIBRARIES ${PBBAM_LIBRARY} ${HTS_LIBRARY})

# set this compiler definition to prevent unhygenic seqan behavior
add_definitions(-DSEQAN_NO_GLOBAL_EXCEPTION_HANDLER)

# shared CXX flags for src & tests
include(CheckCXXCompilerFlag)
set(PacBioCCS_CXX_FLAGS "-std=c++11 -Wall -pedantic -msse3")

# silence gmock warning
check_cxx_compiler_flag("-Wno-unused-variable" HAS_NO_UNUSED_VARIABLE)
if(HAS_NO_UNUSED_VARIABLE)
    set(PacBioCCS_CXX_FLAGS "${PacBioCCS_CXX_FLAGS} -Wno-unused-variable")
endif(HAS_NO_UNUSED_VARIABLE)

# main exe src
add_subdirectory(${ConsensusCore_RootDir} external/ConsensusCore/build)
add_subdirectory(${PBBam_RootDir}         external/pbbam/build)
add_subdirectory(src)

# testing
add_subdirectory(${PacBioCCS_TestsDir})

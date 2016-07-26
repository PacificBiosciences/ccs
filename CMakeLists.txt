##############################################
# CMake build script for the UNANIMITY library
##############################################

cmake_policy(SET CMP0048 NEW)
project(UNANIMITY VERSION 2.0.1 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.2)

set(ROOT_PROJECT_NAME ${PROJECT_NAME} CACHE STRING "root project name")

# Build type
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile RelWithDebInfo" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Build-time options
option(UNY_build_bin    "Build binaries." ON)
option(UNY_build_tests  "Build UNANIMITY's unit tests." ON)
option(UNY_inc_coverage "Include UNANIMITY's coverage script." OFF)
option(UNY_use_ccache   "Build UNANIMITY using ccache, if available." ON)

# Main project paths
set(UNY_RootDir       ${UNANIMITY_SOURCE_DIR})
set(UNY_IncludeDir    ${UNY_RootDir}/include)
set(UNY_SourceDir     ${UNY_RootDir}/src)
set(UNY_SwigDir       ${UNY_RootDir}/swig)
set(UNY_TestsDir      ${UNY_RootDir}/tests)
set(UNY_ThirdPartyDir ${UNY_RootDir}/third-party)

# Project configuration
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})
include(uny-ccache)
include(uny-compilerflags)
include(uny-dependencies)
include(uny-gitsha1)
include(uny-config)

# shared CXX flags for src & tests
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${UNY_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${UNY_DEBUG_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${UNY_LINKER_FLAGS}")

if(APPLE)
    set(PRE_LINK -force_load)
elseif(UNIX)
    SET(PRE_LINK -Wl,-whole-archive)
    SET(POST_LINK -Wl,-no-whole-archive)
endif()

# Build library
add_subdirectory(${UNY_SourceDir})

# Build tests
if(UNY_build_tests)
    add_subdirectory(${UNY_TestsDir})
endif()

# Swig
if (PYTHON_SWIG)
    add_subdirectory(${UNY_SwigDir})
endif()

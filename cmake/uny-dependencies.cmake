
# External libraries

# Get static libraries
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Boost
if(NOT Boost_INCLUDE_DIRS)
    find_package(Boost REQUIRED)
endif()

# pbcopper
if (NOT pbcopper_INCLUDE_DIRS OR
    NOT pbcopper_LIBRARIES)
    if (PYTHON_SWIG)
        set(pbcopper_build_shared OFF CACHE INTERNAL "" FORCE)
    endif()
    set(pbcopper_build_tests OFF CACHE INTERNAL "" FORCE)
    set(pbcopper_build_docs OFF CACHE INTERNAL "" FORCE)
    set(pbcopper_build_examples OFF CACHE INTERNAL "" FORCE)
    add_subdirectory(${UNY_ThirdPartyDir}/pbcopper external/pbcopper/build)
endif()

# only require if NOT called from pip install
if (NOT PYTHON_SWIG)
    # Threads
    if (NOT Threads)
        find_package(Threads REQUIRED)
    endif()

    # ZLIB
    if (NOT ZLIB_INCLUDE_DIRS OR NOT ZLIB_LIBRARIES)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(ZLIB zlib)
    else()
        set(ZLIB_LDFLAGS ${ZLIB_LIBRARIES})
    endif()

    # pbbam
    if (NOT PacBioBAM_INCLUDE_DIRS OR
        NOT PacBioBAM_LIBRARIES)
        set(PacBioBAM_build_docs    OFF CACHE INTERNAL "" FORCE)
        set(PacBioBAM_build_tests   OFF CACHE INTERNAL "" FORCE)
        set(PacBioBAM_build_tools   OFF CACHE INTERNAL "" FORCE)
        add_subdirectory(${UNY_ThirdPartyDir}/pbbam external/pbbam/build)
    endif()

    # cpp-optparse sources
    if (NOT CPPOPTPARSE_CPP)
        set(CPPOPTPARSE_CPP ${UNY_ThirdPartyDir}/cpp-optparse/OptionParser.cpp CACHE INTERNAL "" FORCE)
    endif()

    if (NOT CPPOPTPARSE_IncludeDir)
        set(CPPOPTPARSE_IncludeDir ${UNY_ThirdPartyDir}/cpp-optparse CACHE INTERNAL "" FORCE)
    endif()

    # seqan headers
    if (NOT SEQAN_INCLUDE_DIRS)
        set(SEQAN_INCLUDE_DIRS ${UNY_ThirdPartyDir}/seqan/include CACHE INTERNAL "" FORCE)
    endif()

    # Complete-Striped-Smith-Waterman-Library
    set(ssw_INCLUDE_DIRS ${UNY_ThirdPartyDir}/cssw)
endif()

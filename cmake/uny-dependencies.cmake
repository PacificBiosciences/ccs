
# External libraries

# Get static libraries
SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Boost
if(NOT Boost_INCLUDE_DIRS)
    find_package(Boost REQUIRED)
endif()

# only require if NOT called from pip install
if (NOT PYTHON_SWIG)
    # Threads
    if (NOT Threads)
        find_package(Threads REQUIRED)
    endif()

    # ZLIB
    if (NOT ZLIB_INCLUDE_DIRS OR NOT ZLIB_LIBRARIES)
        find_package(ZLIB REQUIRED)
    endif()

    # htslib
    if(NOT HTSLIB_INCLUDE_DIRS OR
       NOT HTSLIB_LIBRARIES)
        find_package(htslib
                     PATHS ${UNY_ThirdPartyDir}/htslib
                     REQUIRED)
    endif()

    # pbbam
    if (NOT PacBioBAM_INCLUDE_DIRS OR
        NOT PacBioBAM_LIBRARIES)
        set(PacBioBAM_build_docs    OFF CACHE INTERNAL "" FORCE)
        set(PacBioBAM_build_tests   OFF CACHE INTERNAL "" FORCE)
        set(PacBioBAM_build_tools   OFF CACHE INTERNAL "" FORCE)
        add_subdirectory(${UNY_ThirdPartyDir}/pbbam external/pbbam/build)
        if (TARGET htslib)
            add_dependencies(pbbam htslib)
        endif()
    endif()

    # seqan headers
    if (NOT SEQAN_INCLUDE_DIRS)
        set(SEQAN_INCLUDE_DIRS ${UNY_ThirdPartyDir}/seqan/include CACHE INTERNAL "" FORCE)
    endif()

    # pbcopper
    if (NOT pbcopper_INCLUDE_DIRS OR 
        NOT pbcopper_LIBRARIES)
        set(pbcopper_build_tests OFF CACHE INTERNAL "" FORCE)
        set(pbcopper_build_docs OFF CACHE INTERNAL "" FORCE)
        set(pbcopper_build_examples OFF CACHE INTERNAL "" FORCE)
        add_subdirectory(${UNY_ThirdPartyDir}/pbcopper external/pbcopper/build)
    endif()
endif()

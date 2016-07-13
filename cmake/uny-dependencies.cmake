
# External libraries

# pbcopper
if (NOT pbcopper)
	add_subdirectory(${UNY_ThirdPartyDir}/pbcopper external/pbcopper/build)
endif()

# Threads
if (NOT Threads)
	find_package(Threads REQUIRED)
endif()

# Boost
if(NOT Boost_INCLUDE_DIRS)
    find_package(Boost REQUIRED)
endif()

# ZLIB
if (NOT ZLIB_INCLUDE_DIRS OR NOT ZLIB_LIBRARIES)
    find_package(ZLIB REQUIRED)
endif()

# htslib
if(NOT HTSLIB_INCLUDE_DIRS OR
   NOT HTSLIB_LIBRARIES)
    if(HTSLIB_ROOTDIR)
        find_package(htslib
                     PATHS ${HTSLIB_ROOTDIR}
                     REQUIRED)
    else()
        find_package(htslib
                     PATHS ${UNY_ThirdPartyDir}/htslib
                     REQUIRED)
    endif()
endif()

# pbbam
if (NOT pbbam)
    add_subdirectory(${UNY_ThirdPartyDir}/pbbam external/pbbam/build)
endif()

add_dependencies(pbbam htslib)
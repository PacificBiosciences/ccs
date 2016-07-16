
# External libraries

SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# pbcopper
# if (NOT pbcopper)
# 	add_subdirectory(${UNY_ThirdPartyDir}/pbcopper external/pbcopper/build)
# endif()

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
    find_package(htslib
                 PATHS ${UNY_ThirdPartyDir}/htslib
                 REQUIRED)
endif()

# pbbam
if (NOT pbbam)
    add_subdirectory(${UNY_ThirdPartyDir}/pbbam external/pbbam/build)
    if (TARGET htslib)
        add_dependencies(pbbam htslib)
    endif()
endif()


# cpp-optparse sources
set(CPPOPTPARSE_CPP ${UNY_ThirdPartyDir}/cpp-optparse/OptionParser.cpp)
set(CPPOPTPARSE_H   ${UNY_ThirdPartyDir}/cpp-optparse)
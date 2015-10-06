
include_directories(. ${PacBioCCS_IncludeDir})
include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${PBCONSENSUS_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${PBBAM_INCLUDE_DIRS}
    ${CPPOPTPARSE_IncludeDir}
    ${SEQAN_INCLUDE_DIRS}
)

# flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioCCS_CXX_FLAGS}")

# gperftools support
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,-no_pie")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

# libpbccs.a
add_library(pbccs STATIC ${PacBioCCS_CPP} ${PacBioCCS_H})
set_target_properties(pbccs PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PacBioCCS_LibDir}
)

# ccs
set(CCS_Sources
    main/ccs.cpp
    ${CPPOPTPARSE_CPP}
    ${CPPOPTPARSE_H}
)

add_executable(ccs ${CCS_Sources})
set_target_properties(ccs PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioCCS_BinDir}
)
target_link_libraries(ccs
    pbccs
    ${PBCONSENSUS_LIBRARIES}
    ${PBBAM_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)

# gperftools support
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     target_link_libraries(ccs tcmalloc)
# endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

if (TARGET htslib)
    add_dependencies(ccs htslib)
endif()

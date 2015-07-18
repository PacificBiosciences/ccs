
include_directories(.)

include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${ConsensusCore_IncludeDir}
    ${HTS_IncludeDir}
    ${PBBam_IncludeDir}
    ${CPPOPTPARSE_IncludeDir}
    ${SEQAN_IncludeDir}
    ${PacBioCCS_IncludeDir}
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
)

add_executable(ccs ${CCS_Sources})
set_target_properties(ccs PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioCCS_BinDir}
)
add_dependencies(ccs ConsensusCore)
target_link_libraries(ccs
    pbccs
    ${CONSENSUSCORE_LIBRARIES}
    ${PBBAM_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)

# gperftools support
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_libraries(ccs tcmalloc)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

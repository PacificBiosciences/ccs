
# gmock/gtest
set(GMOCK_RootDir    ${ConsensusCore_RootDir}/external/gmock-1.7.0)
set(GMOCK_IncludeDir ${GMOCK_RootDir})

file(GLOB GMOCK_CC   "${GMOCK_RootDir}/*.cc")

# output directory
file(MAKE_DIRECTORY  ${PacBioCCS_TestsDir}/bin)

include_directories(
    ${PacBioCCS_IncludeDir}
)

# Generate paths for test data
configure_file(
    ${PacBioCCS_TestsDir}/TestData.h.in
    ${PacBioCCS_TestsDir}/TestData.h
)

include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${ConsensusCore_IncludeDir}
    ${HTSLIB_INCLUDE_DIRS}
    ${PacBioBAM_IncludeDir}
    ${CPPOPTPARSE_IncludeDir}
    ${SEQAN_IncludeDir}
    ${GMOCK_IncludeDir}
)

file(GLOB Test_CPP "*.cpp")

set(SOURCES
    ${CPPOPTPARSE_CPP}
    ${GMOCK_CC}
    ${Test_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioCCS_CXX_FLAGS}")

add_executable(test_pbccs EXCLUDE_FROM_ALL ${SOURCES})
set_target_properties(test_pbccs PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioCCS_TestsDir}/bin
)
target_link_libraries(test_pbccs
    pbccs
    ConsensusCore
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)
add_custom_target(check test_pbccs)

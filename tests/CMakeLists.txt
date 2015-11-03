
# gmock/gtest
set(GMOCK_RootDir    ${PacBioConsensus_RootDir}/third-party/gmock-1.7.0)
set(GMOCK_IncludeDir ${GMOCK_RootDir})

file(GLOB           GMOCK_CC    "${GMOCK_RootDir}/*.cc")
file(GLOB_RECURSE   GMOCK_H     "${GMOCK_RootDir}/*.h")

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
    ${PBCONSENSUS_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${PBBAM_INCLUDE_DIRS}
    ${CPPOPTPARSE_IncludeDir}
    ${SEQAN_INCLUDE_DIRS}
    ${GMOCK_IncludeDir}
)

file(GLOB PacBioCCSTests_CPP "*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioCCS_CXX_FLAGS}")

add_executable(test_pbccs EXCLUDE_FROM_ALL
    ${PacBioCCSTests_CPP}
    ${CPPOPTPARSE_CPP}
    ${CPPOPTPARSE_H}
    ${GMOCK_CC}
    ${GMOCK_H}
)
set_target_properties(test_pbccs PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioCCS_TestsDir}/bin
)
target_link_libraries(test_pbccs
    pbccs
    ${PBCONSENSUS_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)

if (${ROOT_PROJECT_NAME} STREQUAL "PacBioCCS")
    add_custom_target(check test_pbccs)
endif()

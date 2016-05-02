
# pthread
find_package(Threads)

# gmock/gtest
set(GMOCK_RootDir    ${PacBioConsensus_RootDir}/third-party/gmock-1.7.0)
set(GMOCK_IncludeDir ${GMOCK_RootDir})

file(GLOB           GMOCK_CC    "${GMOCK_RootDir}/*.cc")
file(GLOB_RECURSE   GMOCK_H     "${GMOCK_RootDir}/*.h")

# output directory
file(MAKE_DIRECTORY  ${PacBioConsensus_TestsDir}/bin)

include_directories(
    ${PacBioConsensus_IncludeDir}
)

# Generate paths for test data
#configure_file(
#    ${PacBioConsensus_TestsDir}/TestData.h.in
#    ${PacBioConsensus_TestsDir}/TestData.h
#)

include_directories(SYSTEM
    ${Boost_INCLUDE_DIRS}
    ${GMOCK_IncludeDir}
)

file(GLOB PacBioConsensusTests_CPP "*.cpp")

add_executable(test_pbconsensus EXCLUDE_FROM_ALL
    ${PacBioConsensusTests_CPP}
    ${GMOCK_CC}
    ${GMOCK_H}
)
set_target_properties(test_pbconsensus PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioConsensus_TestsDir}/bin
)
target_link_libraries(test_pbconsensus
    ${PBCONSENSUS_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)

if (${ROOT_PROJECT_NAME} STREQUAL "PacBioConsensus")
    add_custom_target(check test_pbconsensus)
endif()

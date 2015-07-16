
set(GMOCK_RootDir    ${PacBioCCS_ThirdPartyDir}/gmock-1.7.0)
set(GMOCK_IncludeDir ${GMOCK_RootDir})

file(MAKE_DIRECTORY  ${PacBioCCS_TestsDir}/bin)

file(GLOB GMOCK_CC   "${GMOCK_RootDir}/*.cc")

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
    ${HTS_IncludeDir}
    ${PBBam_IncludeDir}
    ${CPPOPTPARSE_IncludeDir}
    ${SEQAN_IncludeDir}
    ${GMOCK_IncludeDir}
)

file(GLOB Test_CPP "*.cpp")

set(SOURCES
    ${CPPOPTPARSE_CPP}
    ${GMOCK_CC}
    ${PacBioCCS_CPP}
    ${Test_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioCCS_CXX_FLAGS}")

add_executable(tests-runner EXCLUDE_FROM_ALL ${SOURCES})
set_target_properties(tests-runner PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioCCS_TestsDir}/bin
)
target_link_libraries(tests-runner
    ${CONSENSUSCORE_LIBRARIES}
    ${PBBAM_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)
add_custom_target(check tests-runner)

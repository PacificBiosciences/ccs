enable_testing()

# pthread
find_package(Threads)

if (EXTENSIVE_TESTING)
    add_definitions(-DEXTENSIVE_TESTING=1)
endif()

# gmock/gtest
set(GMOCK_RootDir    ${UNY_ThirdPartyDir}/gmock-1.7.0)
set(GMOCK_IncludeDir ${GMOCK_RootDir})

file(GLOB           GMOCK_CC    "${GMOCK_RootDir}/*.cc")
file(GLOB_RECURSE   GMOCK_H     "${GMOCK_RootDir}/*.h")

# Generate paths for test data
configure_file(
   ${UNY_TestsDir}/unit/TestData.h.in
   ${CMAKE_BINARY_DIR}/generated/TestData.h
)

include_directories(SYSTEM
    ${GMOCK_IncludeDir}
    ${UNY_ThirdPartyDir}
    ${UNY_IncludeDir}
    ${CMAKE_BINARY_DIR}/generated
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${UNY_ThirdPartyDir}/seqan/include
    unit/
    ${PacBioBAM_INCLUDE_DIRS}
)

file(GLOB UNY_TEST_CPP "unit/*.cpp")

add_executable(test_unanimity EXCLUDE_FROM_ALL
    ${UNY_TEST_CPP}
    ${GMOCK_CC}
    ${GMOCK_H}
)

target_link_libraries(test_unanimity
    ${UNY_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
    ${ZLIB_LIBRARIES}
)

if(${ROOT_PROJECT_NAME} STREQUAL "UNANIMITY")
    add_custom_target(check
        COMMAND ${UNY_RootDir}/tools/check-formatting --all
        COMMAND test_unanimity --gtest_output=xml:${CMAKE_BINARY_DIR}/uny-unit.xml
        COMMAND __PBTEST_CCS_EXE=${CMAKE_BINARY_DIR}/ccs ${UNY_RootDir}/scripts/cram --xunit-file=${CMAKE_BINARY_DIR}/uny-cram.xml --verbose ${UNY_TestsDir}/cram/*.t
        COMMAND __PBTEST_CCS_EXE=${CMAKE_BINARY_DIR}/ccs ${UNY_TestsDir}/python/test_tool_contract.py
        DEPENDS ccs
        WORKING_DIRECTORY ${UNY_RootDir})
endif()

# add_test(test test_unanimity)

if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND UNY_inc_coverage)
    include(uny-codecoverage)
    SETUP_TARGET_FOR_COVERAGE_COBERTURA(coverage test_unanimity coverage)
endif()

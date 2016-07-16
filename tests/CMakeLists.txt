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
    ${CPPOPTPARSE_H}
    unit/
    ${PacBioBAM_INCLUDE_DIRS}
)

file(GLOB UNY_TEST_CPP "unit/*.cpp")

add_executable(test_unanimity
    ${UNY_TEST_CPP}
    ${GMOCK_CC}
    ${GMOCK_H}
)

target_link_libraries(test_unanimity
    ${UNANIMITY_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)

if(${ROOT_PROJECT_NAME} STREQUAL "UNANIMITY")
    add_custom_target(checky test_unanimity
        COMMAND __PBTEST_CCS_EXE=${CMAKE_BINARY_DIR}/ccs ${UNY_RootDir}/scripts/cram.py ${UNY_TestsDir}/cram/*.t
        COMMAND __PBTEST_CCS_EXE=${CMAKE_BINARY_DIR}/ccs ${UNY_TestsDir}/scripts/test_tool_contract.py
        DEPENDS ccs)
endif()

# add_test(test test_unanimity)

file(GLOB           GMOCK_CC    "${GMOCK_RootDir}/*.cc")
file(GLOB_RECURSE   GMOCK_H     "${GMOCK_RootDir}/*.h")

# output directory
file(MAKE_DIRECTORY  ${PacBioChimera_TestsDir}/bin)

include_directories(
    ${PacBioChimera_IncludeDir}
)

include_directories(SYSTEM
    ${GMOCK_IncludeDir}
)

file(GLOB PacBioChimeraTests_CPP "*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioChimera_CXX_FLAGS}")

add_executable(test_pbchimera EXCLUDE_FROM_ALL
    ${PacBioChimeraTests_CPP}
    ${GMOCK_CC}
    ${GMOCK_H}
)

set_target_properties(test_pbchimera PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioChimera_TestsDir}/bin
)

target_link_libraries(test_pbchimera
    ${CMAKE_THREAD_LIBS_INIT}
    ${CMAKE_DL_LIBS}
)

if (${ROOT_PROJECT_NAME} STREQUAL "PacBioChimera")
    add_custom_target(check test_pbchimera)
endif()

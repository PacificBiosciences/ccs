# Config generation

find_git_sha1(UNANIMITY_GIT_SHA1)

file (STRINGS "${UNY_RootDir}/CHANGELOG.md" UNANIMITY_CHANGELOG)

configure_file(
    ${UNY_SourceDir}/UnanimityGitHash.cpp.in
    ${CMAKE_BINARY_DIR}/generated/UnanimityGitHash.cpp
)
configure_file(
    ${UNY_SourceDir}/UnanimityVersion.cpp.in
    ${CMAKE_BINARY_DIR}/generated/UnanimityVersion.cpp
)

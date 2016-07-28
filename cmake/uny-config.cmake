# Config generation

find_git_sha1(UNANIMITY_GIT_SHA1)

file (STRINGS "${UNY_RootDir}/CHANGELOG.md" UNANIMITY_CHANGELOG)

configure_file(
    ${UNY_IncludeDir}/pacbio/Version.h.in
    ${CMAKE_BINARY_DIR}/generated/pacbio/Version.h
)

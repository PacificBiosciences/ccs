# Config generation

find_git_sha1(UNY_GIT_SHA1)

file (STRINGS "${UNY_RootDir}/CHANGELOG.md" UNANIMITY_CHANGELOG)
file (STRINGS "${UNY_RootDir}/version/ccs" CCS_VERSION)

configure_file(
    ${UNY_IncludeDir}/pacbio/Version.h.in
    ${CMAKE_BINARY_DIR}/generated/pacbio/Version.h
)
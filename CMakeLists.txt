########################################################################
# CMake build script for Sparse Alignment library.
########################################################################

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 OLD)
endif()

project(PacBioSparse CXX C)
cmake_minimum_required(VERSION 2.8)

# project version
set(PacBioSparse_VERSION_MAJOR 0)
set(PacBioSparse_VERSION_MINOR 0)
set(PacBioSparse_VERSION_PATCH 1)
set(PacBioSparse_VERSION
  "${PacBioSparse_VERSION_MAJOR}.${PacBioSparse_VERSION_MINOR}.${PacBioSparse_VERSION_PATCH}"
)

# Main project paths
set(PacBioSparse_RootDir        ${PacBioSparse_SOURCE_DIR})
set(PacBioSparse_IncludeDir     ${PacBioSparse_RootDir}/include)

# ccs
set(PacBioCCS_RootDir           ${PacBioSparse_RootDir}/../pbccs CACHE PATH "pbccs root path")

# seqan
set(SeqAn_IncludeDir            ${PacBioSparse_RootDir}/../seqan/include)
link_directories(SeqAn_IncludeDir)

# Declare the PbSparse library target
if (NOT TARGET SPARSE_ALIGNMENT_LIB)
    add_library(SPARSE_ALIGNMENT_LIB INTERFACE)

    target_include_directories(SPARSE_ALIGNMENT_LIB
        INTERFACE
        ${SeqAn_IncludeDir}
        ${PacBioSparse_IncludeDir}
    )
endif()

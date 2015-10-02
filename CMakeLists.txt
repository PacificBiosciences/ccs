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

# SeqAn
if (NOT SeqAn_RootDir AND NOT SeqAn_IncludeDir)
    # Declare the two possible paths to SeqAn
    set(Local_SeqAn_Root        ${PacBioSparse_RootDir}/../seqan)
    set(Biofx_SeqAn_Root        ${PacBioSparse_RootDir}/../../../third-party/cpp/seqan-library-2.0.0)

    # Pick the one that exists or quit
    if (EXISTS ${Local_SeqAn_Root}     AND IS_DIRECTORY ${Local_SeqAn_Root})
        set(SeqAn_RootDir       ${Local_SeqAn_Root})
    elseif (EXISTS ${Biofx_SeqAn_Root} AND IS_DIRECTORY ${Biofx_SeqAn_Root})
        set(SeqAn_RootDir       ${Biofx_SeqAn_Root})
    else (EXISTS ${Local_SeqAn_Root}   AND IS_DIRECTORY ${Local_SeqAn_Root})
        message(FATAL_ERROR     "ERROR: No valid SeqAn path found, exiting")
    endif (EXISTS ${Local_SeqAn_Root}  AND IS_DIRECTORY ${Local_SeqAn_Root})
endif (NOT SeqAn_RootDir AND NOT SeqAn_IncludeDir)

if (NOT SeqAn_IncludeDir)
    set(SeqAn_IncludeDir        ${SeqAn_RootDir}/include)
endif(NOT SeqAn_IncludeDir)

# Declare the PbSparse library target
if (NOT TARGET SPARSE_ALIGNMENT_LIB)
    add_library(SPARSE_ALIGNMENT_LIB INTERFACE)

    target_include_directories(SPARSE_ALIGNMENT_LIB
        INTERFACE
        ${SeqAn_IncludeDir}
        ${PacBioSparse_IncludeDir}
    )
endif()

########################################################################
# CMake build script for Sparse Alignment library.
########################################################################

project(PacBioSparse CXX C)
cmake_minimum_required(VERSION 3.2)

# project version
set(PacBioSparse_MAJOR_VERSION 0)
set(PacBioSparse_MINOR_VERSION 0)
set(PacBioSparse_PATCH_VERSION 1)
set(PacBioSparse_VERSION
  "${PacBioSparse_MAJOR_VERSION}.${PacBioSparse_MINOR_VERSION}.${PacBioSparse_PATCH_VERSION}"
)

# Main project paths
set(PacBioSparse_RootDir        ${PacBioSparse_SOURCE_DIR})
set(PacBioSparse_IncludeDir     ${PacBioSparse_RootDir}/include)

# ccs
set(PacBioCCS_RootDir           ${PacBioLAA_RootDir}/../pbccs                                      CACHE PATH "pbccs root path")
set(PacBioCCS_ThirdPartyDir     ${PacBioCCS_RootDir}/third-party)

# third-party libraries (from pbccs)
set(SeqAn_IncludeDir            ${PacBioCCS_ThirdPartyDir}/seqan-library-2.0.0/include)

# Third-party paths
link_directories(SeqAn_IncludeDir)

if (NOT TARGET PBSPARSE_LIB)
    add_library(PBSPARSE_LIB INTERFACE)

    target_include_directories(PBSPARSE_LIB
		    INTERFACE
        ${SeqAn_IncludeDir}
        ${PacBioSparse_IncludeDir}
    )
endif()

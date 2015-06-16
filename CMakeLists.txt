########################################################################
# CMake build script for Sparse Alignment library.
########################################################################

project(SparseAlignment CXX C)
cmake_minimum_required(VERSION 3.2)

# project version
set(SparseAlignment_MAJOR_VERSION 0)
set(SparseAlignment_MINOR_VERSION 0)
set(SparseAlignment_PATCH_VERSION 1)
set(SparseAlignment_VERSION
  "${SparseAlignment_MAJOR_VERSION}.${SparseAlignment_MINOR_VERSION}.${SparseAlignment_PATCH_VERSION}"
)

# Main project paths
set(SparseAlignment_RootDir       ${SparseAlignment_SOURCE_DIR})
set(SparseAlignment_IncludeDir    ${SparseAlignment_RootDir}/include)
set(Bioinformatics_ThirdPartyDir  ${SparseAlignment_RootDir}/../../third-party/cpp)

# Third-party paths
set(SeqAn_IncludeDir              ${Bioinformatics_ThirdPartyDir}/seqan-library-2.0.0/include)
link_directories(SeqAn_IncludeDir)

if (NOT TARGET SPARSE_ALIGNMENT_LIB)
    add_library(SPARSE_ALIGNMENT_LIB INTERFACE)

    target_include_directories(SPARSE_ALIGNMENT_LIB
		INTERFACE
        ${SeqAn_IncludeDir}
		${SparseAlignment_IncludeDir}
    )
endif()

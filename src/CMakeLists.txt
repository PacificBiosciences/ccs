# get all header files for IDE support
file(GLOB_RECURSE UNY_HEADER "${UNY_IncludeDir}/*.h")
file(GLOB_RECURSE UNY_HIDDEN_HEADER "*.h")
# get sources for src/align
file(GLOB UNY_ALIGN_CPP  "align/*.cpp")
# get sources for src/matrix
file(GLOB UNY_MATRIX_CPP "matrix/*.cpp")
# get sources for src/models
file(GLOB UNY_MODELS_CPP "models/*.cpp")
# get sources for src/poa
file(GLOB UNY_POA_CPP    "poa/*.cpp")

# includes
set(UNY_INCLUDE_DIRS
    ${UNY_IncludeDir}
    ${CMAKE_BINARY_DIR}/generated
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${SEQAN_INCLUDE_DIRS}
    ${CPPOPTPARSE_IncludeDir}
    ${PacBioBAM_INCLUDE_DIRS}
    ${pbcopper_INCLUDE_DIRS}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

# add subset of sources that only contains cc2
add_library(cc2 STATIC
    ${UNY_ALIGN_CPP}
    ${UNY_MATRIX_CPP}
    ${UNY_MODELS_CPP}
    ${UNY_POA_CPP}
    AbstractIntegrator.cpp
    Coverage.cpp
    EasyReadScorer.cpp
    Evaluator.cpp
    EvaluatorImpl.cpp
    ModelConfig.cpp
    ModelFactory.cpp
    ModelFormFactory.cpp
    ModelSelection.cpp
    MonoMolecularIntegrator.cpp
    MultiMolecularIntegrator.cpp
    Mutation.cpp
    Polish.cpp
    Read.cpp
    Recursor.cpp
    Sequence.cpp
    Template.cpp
)

target_include_directories(cc2
    PUBLIC ${UNY_INCLUDE_DIRS}
)

# add main library including everything
add_library(unanimity STATIC
    ${UNY_HEADER}
    ${UNY_HIDDEN_HEADER}
    ${CPPOPTPARSE_CPP}
    AminoAcidCaller.cpp
    ArrayRead.cpp
    ChainSeeds.cpp
    ChemistryMapping.cpp
    ChemistryTriple.cpp
    ConsensusSettings.cpp
    ErrorEstimates.cpp
    Fisher.cpp
    Interval.cpp
    JulietSettings.cpp
    JulietWorkflow.cpp
    MSA.cpp
    MSAColumn.cpp
    ReadId.cpp
    ResistanceCaller.cpp
    SparsePoa.cpp
    SubreadResultCounter.cpp
    Tests.cpp
    Timer.cpp
    Utility.cpp
)

target_include_directories(unanimity
    PUBLIC ${UNY_INCLUDE_DIRS}
)

set(UNY_LIBRARIES
    ${PRE_LINK}
    cc2
    ${POST_LINK}
    unanimity
    ${ZLIB_LIBRARIES}
    ${HTSLIB_LIBRARIES}
    ${PacBioBAM_LIBRARIES}
    ${pbcopper_LIBRARIES}
    CACHE INTERNAL
    "${PROJECT_NAME}: Libraries"
    FORCE
)

# add executables
if (UNY_build_bin)
    ## CCS
    add_executable(ccs ${UNY_SourceDir}/main/ccs.cpp)

    set_target_properties(ccs PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    target_include_directories(ccs
        PUBLIC ${UNY_INCLUDE_DIRS}
    )

    target_link_libraries(ccs
        ${UNY_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        ${CMAKE_DL_LIBS}
    )

    ## Juliet
    add_executable(juliet ${UNY_SourceDir}/main/juliet.cpp)

    set_target_properties(juliet PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    target_include_directories(juliet
        PUBLIC ${UNY_INCLUDE_DIRS}
        PUBLIC ${JUNCTION_ALL_INCLUDE_DIRS}
    )

    target_link_libraries(juliet
        ${UNY_LIBRARIES}
        ${JUNCTION_ALL_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
        ${CMAKE_DL_LIBS}
    )

    if (TARGET pbbam)
        add_dependencies(ccs pbbam)
        add_dependencies(juliet pbbam)
    endif()

    install(TARGETS ccs juliet RUNTIME DESTINATION bin)
endif()

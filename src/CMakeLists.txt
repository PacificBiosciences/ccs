# functions
function(create_exe libName)
    add_executable(${libName} ${UNY_SourceDir}/main/${libName}.cpp)

    set_target_properties(${libName} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    )

    target_link_libraries(${libName}
        ${CMAKE_THREAD_LIBS_INIT}
        ${CMAKE_DL_LIBS}
        ${UNY_LIBRARIES}
    )

    install(TARGETS ${libName} RUNTIME DESTINATION bin)
endfunction()

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
    ${ssw_INCLUDE_DIRS}
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
set(UNANIMITY_FILES
    ${UNY_HEADER}
    ${UNY_HIDDEN_HEADER}
    ${CPPOPTPARSE_CPP}
    ChemistryMapping.cpp
    ChemistryTriple.cpp
    Interval.cpp
    ReadId.cpp
    SparsePoa.cpp
    SubreadResultCounter.cpp
    Timer.cpp
    Utility.cpp
)

if (NOT PYTHON_SWIG)
    set(UNANIMITY_FILES ${UNANIMITY_FILES}
        ${ssw_INCLUDE_DIRS}/ssw_cpp.cpp
        ${ssw_INCLUDE_DIRS}/ssw.c
        ConsensusSettings.cpp
        LocalAlignment.cpp
    )
endif()

add_library(unanimity STATIC
    ${UNANIMITY_FILES}
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
    create_exe(ccs)
endif()

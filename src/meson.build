#############
# libraries #
#############

uny_cc2_cpp_sources = [
  # -------
  # align
  # -------
  'align/AffineAlignment.cpp',
  'align/AlignConfig.cpp',
  'align/BandedChainAlignment.cpp',
  'align/LinearAlignment.cpp',
  'align/PairwiseAlignment.cpp',

  # --------
  # matrix
  # --------
  'matrix/BasicDenseMatrix.cpp',
  'matrix/ScaledMatrix.cpp',
  'matrix/SparseMatrix.cpp',

  # --------
  # models
  # --------
  'models/MarginalModelForm.cpp',
  'models/P6C4NoCovModel.cpp',
  'models/PwSnrAModelForm.cpp',
  'models/PwSnrModelForm.cpp',
  'models/S_P1C1Beta_Model.cpp',
  'models/S_P1C1v1_Model.cpp',
  'models/S_P1C1v2_Model.cpp',
  'models/S_P2C2v5_Model.cpp',
  'models/SnrModelForm.cpp',

  # -----
  # poa
  # -----
  'poa/PoaAlignmentMatrix.cpp',
  'poa/PoaConsensus.cpp',
  'poa/PoaGraph.cpp',
  'poa/PoaGraphImpl.cpp',
  'poa/PoaGraphTraversals.cpp',
  'poa/RangeFinder.cpp',

  # -----
  # cc2
  # -----
  'Coverage.cpp',
  'EasyReadScorer.cpp',
  'Evaluator.cpp',
  'EvaluatorImpl.cpp',
  'Integrator.cpp',
  'IntervalMask.cpp',
  'ModelConfig.cpp',
  'ModelFactory.cpp',
  'ModelFormFactory.cpp',
  'ModelSelection.cpp',
  'Mutation.cpp',
  'Polish.cpp',
  'PolishResult.cpp',
  'Read.cpp',
  'Recursor.cpp',
  'Sequence.cpp',
  'Template.cpp',

  # --------------
  # cpp-optparse
  # --------------
  '../third-party/cpp-optparse/OptionParser.cpp',

  # ---------
  # supplib
  # ---------
  'ChemistryMapping.cpp',
  'ChemistryTriple.cpp',
  'Interval.cpp',
  'ReadId.cpp',
  'SparsePoa.cpp',
  'SubreadResultCounter.cpp',
  'Timer.cpp',
  'Utility.cpp',

  'ConsensusSettings.cpp',
  'LocalAlignment.cpp',

  uny_ccsw_cpp_sources]

uny_cc2_lib_impl = library(
  'unanimity',
  uny_cc2_cpp_sources,
  install : true,
  soversion : 0,
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_boost_dep, uny_thread_dep, uny_seqan_dep],
  include_directories : [uny_include_directories, uny_ccsw_header],
  cpp_args : [uny_warning_flags])

if get_option('default_library') == 'shared'
  uny_cc2_lib_shared = uny_cc2_lib_impl
  uny_cc2_lib_static = []
else
  uny_cc2_lib_shared = []
  uny_cc2_lib_static = uny_cc2_lib_impl
endif

###############
# executables #
###############

# ccs
if get_option('enable-uny-build-ccs')
  uny_ccs = executable(
  'ccs', [
    'main/ccs.cpp'],
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_boost_dep, uny_thread_dep],
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib_shared,
  link_whole : uny_cc2_lib_static,
  cpp_args : uny_warning_flags)
endif

# ChimeraLabeler
if get_option('enable-uny-build-chimera')
  executable(
  'ChimeraLabeler', [
    'main/ChimeraLabeler.cpp'],
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_seqan_dep],
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib_shared,
  link_whole : uny_cc2_lib_static,
  cpp_args : uny_warning_flags)
endif

# ccs_sim
if get_option('enable-uny-build-sim')
  executable(
  'ccs_sim', [
    'main/ccs_sim.cpp'],
  dependencies : [uny_pbbam_dep],
  include_directories : [uny_include_directories, uny_ccsw_header],
  link_with : uny_cc2_lib_shared,
  link_whole : uny_cc2_lib_static,
  cpp_args : uny_warning_flags)
endif

###################
# dependency info #
###################

import('pkgconfig').generate(
  libraries : uny_cc2_lib_impl,
  version : meson.project_version(),
  name : 'unanimity',
  filebase : 'unanimity',
  description : 'PacBio C++ library for calling circular consensus sequences')

declare_dependency(
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib_impl,
  version : meson.project_version(),
  sources : uny_cc2_cpp_sources)

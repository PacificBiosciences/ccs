################
# version info #
################
uny_ver_arr = meson.project_version().split('.')
uny_major_version = uny_ver_arr[0]
uny_minor_version = uny_ver_arr[1]
uny_patch_version = uny_ver_arr[2]

uny_config = configuration_data()
uny_config.set('UNANIMITY_VERSION', meson.project_version())
uny_config.set('UNANIMITY_VERSION_MAJOR', uny_major_version)
uny_config.set('UNANIMITY_VERSION_MINOR', uny_minor_version)
uny_config.set('UNANIMITY_VERSION_PATCH', uny_patch_version)

uny_version_sources = [
  vcs_tag(
    input : 'UnanimityGitHash.cpp.in',
    output : 'UnanimityGitHash.cpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@UNANIMITY_GIT_SHA1@'),
  configure_file(
    input : 'UnanimityVersion.cpp.in',
    output : 'UnanimityVersion.cpp',
    configuration : uny_config)]

#############
# libraries #
#############

uny_cc2_cpp_sources = files([
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
  'LocalAlignment.cpp'])

# models
uny_models = [
  # compile-time parameterized
  ['P6C4NoCov', files('models/P6C4NoCov_Model.cpp')],
  ['S_P1C1Beta', files('models/S_P1C1Beta_Model.cpp')],
  ['S_P1C1v1', files('models/S_P1C1v1_Model.cpp')],
  ['S_P1C1v2', files('models/S_P1C1v2_Model.cpp')],
  ['S_P2C2v5', files('models/S_P2C2v5_Model.cpp')],

  # run-time parameterized
  ['Marginal', files('models/MarginalModelForm.cpp')],
  ['PwSnrA', files('models/PwSnrAModelForm.cpp')],
  ['PwSnr', files('models/PwSnrModelForm.cpp')],
  ['Snr', files('models/SnrModelForm.cpp')]]

uny_enabled_models_str = ''
foreach i : uny_models
  uny_enabled_models_str += ' X(' + i[0] + ')'
  uny_cc2_cpp_sources += i[1]
endforeach

uny_enabled_models = configuration_data()
uny_enabled_models.set('ENABLED_MODELS', uny_enabled_models_str)

uny_cc2_cpp_sources += configure_file(
  input : 'ModelInternalInitializer.cpp.in',
  output : 'ModelInternalInitializer.cpp',
  configuration : uny_enabled_models)

uny_cc2_cpp_sources += uny_version_sources
uny_cc2_cpp_sources += uny_cssw_cpp_sources

uny_cc2_lib = library(
  'unanimity',
  uny_cc2_cpp_sources,
  install : true,
  soversion : 0,
  version : meson.project_version(),
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_boost_dep, uny_thread_dep, uny_seqan_dep],
  include_directories : [uny_include_directories, uny_cssw_header],
  cpp_args : uny_warning_flags)

###############
# executables #
###############

# ccs
if get_option('enable-build-ccs')
  uny_ccs = executable(
  'ccs', [
    'main/ccs.cpp'],
  install : true,
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_boost_dep, uny_thread_dep],
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib,
  cpp_args : uny_warning_flags)
endif

# arrow
if get_option('enable-build-arrow')
  uny_arrow = executable(
  'arrow', [
    'main/arrow.cpp'],
  install : true,
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_boost_dep, uny_thread_dep],
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib,
  cpp_args : uny_warning_flags)
endif

# ChimeraLabeler
if get_option('enable-build-chimera')
  executable(
  'ChimeraLabeler', [
    'main/ChimeraLabeler.cpp'],
  install : true,
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_seqan_dep],
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib,
  cpp_args : uny_warning_flags)
endif

# ccs_sim
if get_option('enable-build-sim')
  executable(
  'ccs_sim', [
    'main/ccs_sim.cpp'],
  install : true,
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep],
  include_directories : [uny_include_directories, uny_cssw_header],
  link_with : uny_cc2_lib,
  cpp_args : uny_warning_flags)
endif

# genomic_sim
if get_option('enable-build-sim')
  executable(
  'genomic_sim', [
    'main/genomic_sim.cpp'],
  install : true,
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep],
  include_directories : [uny_include_directories, uny_cssw_header],
  link_with : uny_cc2_lib,
  cpp_args : uny_warning_flags)
endif

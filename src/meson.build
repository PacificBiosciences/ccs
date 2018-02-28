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

  # ------------------
  # genomicconsensus
  # ------------------
  'genomicconsensus/experimental/Consensus.cpp',
  'genomicconsensus/experimental/ConsensusModelFactory.cpp',
  'genomicconsensus/experimental/Filters.cpp',
  'genomicconsensus/experimental/GenomicConsensus.cpp',
  'genomicconsensus/experimental/IPoaModel.cpp',
  'genomicconsensus/experimental/Input.cpp',
  'genomicconsensus/experimental/Intervals.cpp',
  'genomicconsensus/experimental/Output.cpp',
  'genomicconsensus/experimental/Settings.cpp',
  'genomicconsensus/experimental/SettingsOptions.h',
  'genomicconsensus/experimental/SettingsToolContract.h',
  'genomicconsensus/experimental/Sorting.cpp',
  'genomicconsensus/experimental/Workflow.cpp',
  'genomicconsensus/experimental/arrow/ArrowModel.cpp',
  'genomicconsensus/experimental/plurality/PluralityModel.cpp',
  'genomicconsensus/experimental/poa/PoaModel.cpp',

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

  # ------
  # gcpp
  # ------
  'genomicconsensus/experimental/Consensus.cpp',
  'genomicconsensus/experimental/ConsensusModelFactory.cpp',
  'genomicconsensus/experimental/Filters.cpp',
  'genomicconsensus/experimental/GenomicConsensus.cpp',
  'genomicconsensus/experimental/Input.cpp',
  'genomicconsensus/experimental/Intervals.cpp',
  'genomicconsensus/experimental/IPoaModel.cpp',
  'genomicconsensus/experimental/Output.cpp',
  'genomicconsensus/experimental/Settings.cpp',
  'genomicconsensus/experimental/Sorting.cpp',
  'genomicconsensus/experimental/Workflow.cpp',

  'genomicconsensus/experimental/arrow/ArrowModel.cpp',

  'genomicconsensus/experimental/io/FastaWriter.cpp',
  'genomicconsensus/experimental/io/FastqWriter.cpp',
  'genomicconsensus/experimental/io/GffWriter.cpp',
  'genomicconsensus/experimental/io/VcfWriter.cpp',

  'genomicconsensus/experimental/plurality/PluralityModel.cpp',

  'genomicconsensus/experimental/poa/PoaModel.cpp',

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

# install library if
# - either running as a proper project
# - or using shared libraries
uny_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

uny_cc2_lib = library(
  'unanimity',
  uny_cc2_cpp_sources,
  install : uny_lib_install,
  # use boost SONAME practice:
  #   cause ld startup issues before
  #   you even have the chance of running
  #   into ABI issues.
  soversion : meson.project_version(),
  version : meson.project_version(),
  dependencies : [
    uny_pbcopper_dep,
    uny_pbbam_dep,
    uny_boost_dep,
    uny_thread_dep,
    uny_seqan_dep],
  include_directories : [
    uny_include_directories,
    uny_cssw_header],
  cpp_args : uny_warning_flags)

###############
# executables #
###############

if not meson.is_subproject()
  # ccs
  uny_ccs = executable(
    'ccs', files([
      'main/ccs.cpp']),
    install : true,
    dependencies : [
      uny_pbcopper_dep,
      uny_pbbam_dep,
      uny_boost_dep,
      uny_thread_dep,
      uny_zlib_dep],
    include_directories : uny_include_directories,
    link_with : uny_cc2_lib,
    cpp_args : uny_warning_flags)

  # arrow
  uny_arrow = executable(
    'arrow', files([
      'main/arrow.cpp']),
    install : true,
    dependencies : [
      uny_pbcopper_dep,
      uny_pbbam_dep,
      uny_boost_dep,
      uny_thread_dep,
      uny_zlib_dep],
    include_directories : uny_include_directories,
    link_with : uny_cc2_lib,
    cpp_args : uny_warning_flags)

  # ChimeraLabeler
  executable(
    'ChimeraLabeler', files([
      'main/ChimeraLabeler.cpp']),
    install : true,
    dependencies : [
      uny_pbcopper_dep,
      uny_pbbam_dep,
      uny_seqan_dep],
    include_directories : uny_include_directories,
    link_with : uny_cc2_lib,
    cpp_args : uny_warning_flags)

  # ccs_sim
  executable(
    'ccs_sim', files([
      'main/ccs_sim.cpp']),
    install : true,
    dependencies : [
      uny_pbcopper_dep,
      uny_pbbam_dep,
      uny_zlib_dep,
      uny_boost_dep],
    include_directories : [
      uny_include_directories,
      uny_cssw_header],
    link_with : uny_cc2_lib,
    cpp_args : uny_warning_flags)

  # gcpp
  uny_gcpp = executable(
    'gcpp', files([
      'main/gcpp.cpp']),
    install : true,
    dependencies : [
      uny_pbcopper_dep,
      uny_pbbam_dep,
      uny_boost_dep,
      uny_thread_dep,
      uny_zlib_dep],
    include_directories : uny_include_directories,
    link_with : uny_cc2_lib,
    cpp_args : uny_warning_flags)

  # genomic_sim
  executable(
    'genomic_sim', files([
      'main/genomic_sim.cpp']),
    install : true,
    dependencies : [
      uny_pbcopper_dep,
      uny_pbbam_dep,
      uny_zlib_dep,
      uny_boost_dep],
    include_directories : [
      uny_include_directories,
      uny_cssw_header],
    link_with : uny_cc2_lib,
    cpp_args : uny_warning_flags)

  # variantCaller
  uny_variantCaller = executable(
    'variantCaller', files([
      'variantCaller/VariantCallerSettings.cpp',
      'variantCaller/variantCaller.cpp']),
    install : true,
    dependencies : [
      uny_pbcopper_dep,
      uny_pbbam_dep,
      uny_boost_dep,
      uny_thread_dep,
      uny_zlib_dep],
    include_directories : uny_include_directories,
    link_with : uny_cc2_lib,
    cpp_args : uny_warning_flags)
endif

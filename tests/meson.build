uny_test_cpp_sources = [
  'unit/ChimeraTest.cpp',
  'unit/Mutations.cpp',
  'unit/RandomDNA.cpp',
  'unit/TestAlignment.cpp',
  'unit/TestBandedChainAlign.cpp',
  'unit/TestChemistry.cpp',
  'unit/TestConsensus.cpp',
  'unit/TestCoverage.cpp',
  'unit/TestIntegrator.cpp',
  'unit/TestInterval.cpp',
  'unit/TestIntervalMask.cpp',
  'unit/TestLoadModels.cpp',
  'unit/TestMutationEnumerator.cpp',
  'unit/TestPoaConsensus.cpp',
  'unit/TestPolish.cpp',
  'unit/TestSequence.cpp',
  'unit/TestSparseAlign.cpp',
  'unit/TestSparsePoa.cpp',
  'unit/TestSparseVector.cpp',
  'unit/TestTemplate.cpp',
  'unit/TestUtility.cpp',
  'unit/TestWhitelist.cpp'
]

uny_TestData_h_config = configuration_data()
uny_TestData_h_config.set('UNY_TestsDir', meson.current_source_dir())
uny_TestData_h_config.set('CMAKE_BINARY_DIR', meson.current_build_dir())
configure_file(input : 'unit/TestData.h.in', output : 'TestData.h', configuration : uny_TestData_h_config)

# find GoogleTest and GoogleMock
uny_gtest_dep = dependency('gtest', fallback : ['gtest', 'gtest_with_main_dep'])
uny_gmock_dep = dependency('gmock', fallback : ['gtest', 'gmock_mainless_dep'])

uny_test = executable(
  'uny_test',
  [uny_test_cpp_sources],
  dependencies : [
    uny_pbbam_dep,
    uny_pbcopper_dep,
    uny_gtest_dep,
    uny_gmock_dep,
    uny_boost_dep,
    uny_seqan_dep],
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib_shared,
  link_whole : uny_cc2_lib_static,
  cpp_args : [uny_warning_flags],
  install : false)

#########
# tests #
#########

test('unanimity formatting check', uny_clang_formatter, args : ['--all'], workdir : meson.source_root())

test('unanimity testsuite', uny_test, args : ['--gtest_output=xml:uny-unit.xml'])

test(
  'unanimity ccs test',
  uny_cram_script,
  args : [
    '--xunit-file=uny-cram.xml',
    '--verbose'] +
    uny_cram_test_data,
  env : [
    '__PBTEST_CCS_EXE=' + uny_ccs.full_path(),
    'TESTDIR=' + join_paths([meson.current_source_dir(), 'cram'])],
  timeout : 1800)

test(
  'unanimity ccs contract',
  uny_test_tool_contract,
  env : ['__PBTEST_CCS_EXE=' + uny_ccs.full_path()])

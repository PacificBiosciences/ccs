subdir('unit')

uny_TestData_h_config = configuration_data()
uny_TestData_h_config.set('UNY_TestsDir', meson.current_source_dir())
uny_TestData_h_config.set('CMAKE_BINARY_DIR', meson.current_build_dir())
uny_TestData_h = configure_file(
  input : files('unit/TestData.h.in'),
  output : 'TestData.h',
  configuration : uny_TestData_h_config)
uny_test_cpp_sources += uny_TestData_h

# find GoogleTest and GoogleMock
uny_gmock_main_dep = dependency('gmock_main', fallback : ['gtest', 'gmock_main_dep'])

uny_test = executable(
  'uny_test',
  uny_test_cpp_sources,
  dependencies : [
    uny_pbbam_dep,
    uny_pbcopper_dep,
    uny_gmock_main_dep,
    uny_boost_dep,
    uny_seqan_dep],
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib,
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

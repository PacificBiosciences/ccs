# Determine numpy includedir
uny_numpy_incdir = run_command(uny_python_interp, ['-c', 'import sys; import numpy; sys.stdout.write(numpy.get_include())'])
if uny_numpy_incdir.returncode() != 0
  error('Could not determine NumPy includedir!')
else
  message('NumPy includedir: ' + uny_numpy_incdir.stdout())
  uny_numpy_include_directories = include_directories(uny_numpy_incdir.stdout())
endif

#############
# libraries #
#############

uny_consensus_core2_wrap_cpp = custom_target('swig2cpp',
  output : ['ConsensusCore2.py', 'ConsensusCore2_wrap.cpp'],
  input : 'ConsensusCore2.i',
  command : [
    uny_swig,
    '-Wextra',
    '-c++',
    '-python',
    '-builtin',
    '-module',
    'ConsensusCore2',
    '-I' + join_paths([meson.source_root(), 'include']),
    '-I' + join_paths([meson.build_root(), 'include']),
    '-o',
    '@OUTPUT1@',
    '-outdir',
    meson.current_build_dir(),
    '@INPUT@'],
  install : true,
  install_dir : [uny_python_install_dir, false])

# python .py files should be byte-compiled before installing, in order not to have the python
# vm generate them at runtime which on some distros causes them to become orphaned
custom_target('uny_py_bytecompile_pyo',
  input: uny_consensus_core2_wrap_cpp,
  output: 'ConsensusCore2.pyo',
  command: [uny_python_interp, '-OO', '-m', 'compileall', '-q', '-f', '@INPUT0@'],
  install : true,
  install_dir : uny_python_install_dir)

custom_target('uny_py_bytecompile_pyc',
  input: uny_consensus_core2_wrap_cpp,
  output: 'ConsensusCore2.pyc',
  command: [uny_python_interp, '-m', 'compileall', '-q', '-f', '@INPUT0@'],
  install : true,
  install_dir : uny_python_install_dir)

uny_consensus_core2_lib = shared_library(
  '_ConsensusCore2',
  uny_consensus_core2_wrap_cpp,
  name_prefix : '',
  # python is brain-damaged and will
  # not load dylib libraries on macOS
  # https://stackoverflow.com/a/2582504
  name_suffix : 'so',
  install : true,
  install_dir : uny_python_install_dir,
  dependencies : [uny_boost_dep, uny_python_dep],
  include_directories : [
    uny_numpy_include_directories,
    uny_include_directories],
  link_with : uny_cc2_lib,
  cpp_args : uny_swig_warning_flags)

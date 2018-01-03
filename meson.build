project(
  'UNANIMITY',
  ['cpp', 'c'],
  version : '3.0.0',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++14',
    'c_std=c11'],
  license : 'BSD-3',
  meson_version : '>= 0.41.1')

############
# CXXFLAGS #
############

uny_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-unused-parameter',
  '-Wno-non-virtual-dtor',
  '-Wno-unused-variable',
  '-Wno-sign-compare']
    if cpp.has_argument(cflag)
      uny_warning_flags += cflag
    endif
endforeach

if cpp.get_id() == 'clang'
  foreach cflag: [
    '-Wno-unused-local-typedef']
      if cpp.has_argument(cflag)
        uny_warning_flags += cflag
      endif
  endforeach
endif

uny_swig_warning_flags = uny_warning_flags
foreach cflag: [
  '-Wno-delete-non-virtual-dtor',
  '-Wno-dynamic-class-memaccess']
    if cpp.has_argument(cflag)
      uny_swig_warning_flags += cflag
    endif
endforeach

################
# dependencies #
################

# threads
uny_thread_dep = dependency('threads', required : true)

# boost
uny_boost_dep = dependency('boost', required : true)

# pbbam
uny_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])

# pbcopper
uny_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])

# zlib
uny_zlib_dep = dependency('zlib', required : true)

# seqan
uny_seqan_dep = dependency('seqan-2', required : true)

# cssw
uny_cssw_header = include_directories('third-party/cssw/')
uny_cssw_cpp_sources = files([
  'third-party/cssw/ssw_cpp.cpp',
  'third-party/cssw/ssw.c'])

###########
# headers #
###########

subdir('include')

#####################
# sources + library #
#####################

subdir('src')

########
# swig #
########

if (not meson.is_subproject()) and get_option('swig')
  uny_python_interp = find_program(get_option('PYTHON'), required : true)
  uny_swig = find_program('swig', required : true)

  uny_pv_ver = run_command(uny_python_interp, ['-c', 'import sys; sys.stdout.write(sys.version[:3])']).stdout()
  message('Python version: ' + uny_pv_ver)

  if get_option('pythoninstalldir') == ''
    # user did not provide python install dir,
    # retrieving from python interpreter.
    uny_python_sitedir_config = configuration_data()
    uny_python_sitedir_config.set('PREFIX', get_option('prefix'))
    uny_python_sitedir_script = configure_file(
      input : 'scripts/meson/python_sitedir.py.in',
      output : 'python_sitedir.py',
      configuration : uny_python_sitedir_config)

    uny_python_install_dir = run_command(uny_python_interp, uny_python_sitedir_script).stdout()
  else
    # user provided python install dir, do not
    # override it.
    uny_python_install_dir = get_option('pythoninstalldir')
  endif

  message('Python modules will be installed into: ' + uny_python_install_dir)

  # python2 compiling + linking
  uny_python_dep = dependency('python-2.7', required : true)

  subdir('swig')
endif

#########
# tests #
#########

if (not meson.is_subproject()) and get_option('tests')
  uny_clang_formatter = find_program('tools/check-formatting')
  uny_cram_script = find_program('scripts/cram')
  uny_cram_test_data = files([
    'tests/cram/0passes.t',
    'tests/cram/100zmws.t',
    'tests/cram/100zmws_byStrand.t',
    'tests/cram/43385286.t',
    'tests/cram/force.t',
    'tests/cram/logfile.t',
    'tests/cram/nozfilter.t',
    'tests/cram/tiny.t',
    'tests/cram/tiny_noPolish.t'])
  uny_test_tool_contract = find_program('tests/python/test_tool_contract.py')

  subdir('tests')
endif

###################
# dependency info #
###################

if (not meson.is_subproject())
  import('pkgconfig').generate(
    libraries : uny_cc2_lib,
    version : meson.project_version(),
    name : 'unanimity',
    filebase : 'unanimity',
    description : 'PacBio C++ library for calling circular consensus sequences')
endif

uny_dep = declare_dependency(
  include_directories : uny_include_directories,
  link_with : uny_cc2_lib,
  dependencies : [uny_pbcopper_dep, uny_pbbam_dep, uny_boost_dep, uny_thread_dep, uny_seqan_dep],
  version : meson.project_version(),
  compile_args : uny_swig_warning_flags)

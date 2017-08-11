# replace version strings
uny_ver_arr = meson.project_version().split('.')
uny_major_version = uny_ver_arr[0]
uny_minor_version = uny_ver_arr[1]
uny_patch_version = uny_ver_arr[2]

uny_config = configuration_data()
uny_config.set('UNANIMITY_VERSION', meson.project_version())
uny_config.set('UNANIMITY_VERSION_MAJOR', uny_major_version)
uny_config.set('UNANIMITY_VERSION_MINOR', uny_minor_version)
uny_config.set('UNANIMITY_VERSION_PATCH', uny_patch_version)

configure_file(
  input : 'UnanimityVersion.h.in',
  output : 'UnanimityVersion.h',
  configuration : uny_config)

# replace git commit id
vcs_tag(
  input : 'UnanimityGitHash.h.in',
  output : 'UnanimityGitHash.h',
  command : ['git', 'describe', '--always', '--dirty=*'],
  replace_string : '@UNANIMITY_GIT_SHA1@')

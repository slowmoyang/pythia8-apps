project(
  'pythia8-apps',
  'cpp',
  license: 'MIT',
  version: '1.0.0',
)

# root
root_config = find_program('root-config', required: true)
root_compile_args = run_command(root_config, '--cflags', check: true).stdout().strip().split(' ')
root_link_args = run_command(root_config, '--ldflags', '--libs', check: true).stdout().strip().split(' ')
root_include_directories = run_command(root_config, '--incdir', check: true).stdout().strip().split(' ')
root_dep = declare_dependency(
    compile_args: root_compile_args,
    link_args: root_link_args,
    include_directories: root_include_directories,
)

# pythia8
pythia8_config = find_program('pythia8-config', required: true)
pythia8_compile_args = run_command(pythia8_config, '--cxxflags', check: true).stdout().strip().split(' ')
pythia8_link_args = run_command(pythia8_config, '--ldflags', '--libs', check: true).stdout().strip().split(' ')
pythia8_include_directories = run_command(pythia8_config, '--includedir', check: true).stdout().strip().split(' ')
pythia8_dep = declare_dependency(
    compile_args: pythia8_compile_args,
    link_args: pythia8_link_args,
    include_directories: pythia8_include_directories,
)

# hepmc3
hepmc3_config = find_program('HepMC3-config', required: true)
hepmc3_compile_args = run_command(hepmc3_config, '--cppflags', check: true).stdout().strip().split(' ')
hepmc3_link_args = run_command(hepmc3_config, '--ldflags', check: true).stdout().strip().split(' ')
hepmc3_include_directories = run_command(hepmc3_config, '--includedir', check: true).stdout().strip().split(' ')
hepmc3_dep = declare_dependency(
    compile_args: hepmc3_compile_args,
    link_args: hepmc3_link_args,
    include_directories: hepmc3_include_directories,
)

argparse_dep = dependency(
  'argparse',
  required: true,
)

dependencies = [
  root_dep,
  pythia8_dep,
  hepmc3_dep,
  argparse_dep,
]

executable(
  'pythia8-hadronizer',
  './src/pythia8-hadronizer.cc',
  dependencies: dependencies,
  install : true,
)

executable(
  'minbias-generation',
  './src/minbias-generation.cc',
  dependencies: dependencies,
  install: true,
)

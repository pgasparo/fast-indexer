# Meson file for FDIP
project('libxgandalf', ['cpp'],
        version: '1.0',
        meson_version: '>=0.49.0',
        default_options: ['buildtype=debugoptimized'])

if get_option('buildtype') == 'custom'
  add_project_arguments('-O3', '-flto',  language : 'cpp')
  add_project_link_arguments('-O3', '-flto',  language : 'cpp')
endif

extra_cflag_list = []

compiler = meson.get_compiler('cpp')
if get_option('native')
        if compiler.get_id() == 'gcc'
                message('building for native CPU (gcc)')
                add_project_arguments('-march=native', language: 'cpp')
                extra_cflag_list += '-march=native'
        elif compiler.get_id() == 'clang'
                message('building for native CPU (clang)')
                add_project_arguments('-march=native', language: 'cpp')
                extra_cflag_list += '-march=native'
        endif
endif
if get_option('omp')
        if compiler.get_id() == 'gcc'
                message('building with openmp (gcc)')
                add_project_arguments('-fopenmp', language: 'cpp')
                extra_cflag_list += '-fopenmp'
                add_project_link_arguments('-fopenmp', language: 'cpp')
        elif compiler.get_id() == 'clang'
                message('building with openmp (clang)')
                add_project_arguments('-fopenmp', language: 'cpp')
                extra_cflag_list += '-fopenmp'
                add_project_link_arguments('-fopenmp', language: 'cpp')
        endif
endif

incdir = include_directories('include')

if get_option('download_eigen')
        eigen_proj = subproject('eigen')
        eigendep = eigen_proj.get_variable('eigen_dep')
else
        eigendep = dependency('eigen3',
                        required: true,
                        fallback: ['eigen', 'eigen_dep'])
endif

add_project_arguments('-DEIGEN_NO_DEBUG', '-DEIGEN_NO_AUTOMATIC_RESIZING', language: 'cpp')
extra_cflag_list += [ '-DEIGEN_NO_DEBUG', '-DEIGEN_NO_AUTOMATIC_RESIZING' ]

if get_option('mkl')
        mkldep = dependency('mkl-dynamic-lp64-' + get_option('mkl_flavour'), required: false)
        if mkldep.found()
                if get_option('mkl_direct')
                        message('using MKL direct')
                        add_project_arguments('-DEIGEN_USE_MKL', language: 'cpp')
                        extra_cflag_list += '-DEIGEN_USE_MKL'
                else
                        message('using MKL')
                endif
                add_project_arguments('-DEIGEN_USE_MKL_ALL', language: 'cpp')
                extra_cflag_list += '-DEIGEN_USE_MKL_ALL'
        endif
else
        mkldep = dependency('', required: false)
endif

# Required for xgandalf
precomputed_point_dir = get_option('prefix')/get_option('datadir')/'xgandalf/precomputedSamplePoints'
add_project_arguments('-DPRECOMPUTED_DATA_DIR='+precomputed_point_dir, language: 'cpp')

xgandalf = library('xgandalf',
                   ['src/Dbscan.cpp',
                    'src/DetectorToReciprocalSpaceTransform.cpp',
                    'src/ExperimentSettings.cpp',
                    'src/HillClimbingOptimizer.cpp',
                    'src/IndexerAutocorrPrefit.cpp',
                    'src/IndexerBase.cpp',
                    'src/IndexerPlain.cpp',
                    'src/InverseSpaceTransform.cpp',
                    'src/LatticeAssembler.cpp',
                    'src/Lattice.cpp',
                    'src/pointAutocorrelation.cpp',
                    'src/ReciprocalToRealProjection.cpp',
                    'src/refinement.cpp',
                    'src/samplePointsFiltering.cpp',
                    'src/SamplePointsGenerator.cpp',
                    'src/SimpleMonochromaticDiffractionPatternPrediction.cpp',
                    'src/SimpleMonochromaticProjection.cpp',
                    'src/SparsePeakFinder.cpp'],
                   include_directories: [incdir],
                   dependencies: [eigendep, mkldep],
                   install: true,
                   override_options : ['cpp_std=c++17'])

install_headers(['include/xgandalf/BadInputException.h',
                 'include/xgandalf/CustomException.h',
                 'include/xgandalf/Dbscan.h',
                 'include/xgandalf/DetectorToReciprocalSpaceTransform.h',
                 'include/xgandalf/eigenDiskImport.h',
                 'include/xgandalf/eigenSTLContainers.h',
                 'include/xgandalf/ExperimentSettings.h',
                 'include/xgandalf/HillClimbingOptimizer.h',
                 'include/xgandalf/IndexerAutocorrPrefit.h',
                 'include/xgandalf/IndexerBase.h',
                 'include/xgandalf/IndexerPlain.h',
                 'include/xgandalf/InverseSpaceTransform.h',
                 'include/xgandalf/LatticeAssembler.h',
                 'include/xgandalf/Lattice.h',
                 'include/xgandalf/pointAutocorrelation.h',
                 'include/xgandalf/ReciprocalToRealProjection.h',
                 'include/xgandalf/refinement.h',
                 'include/xgandalf/samplePointsFiltering.h',
                 'include/xgandalf/SamplePointsGenerator.h',
                 'include/xgandalf/SimpleMonochromaticDiffractionPatternPrediction.h',
                 'include/xgandalf/SimpleMonochromaticProjection.h',
                 'include/xgandalf/SparsePeakFinder.h',
                 'include/xgandalf/WrongUsageException.h'],
                subdir: 'xgandalf')

install_data(['precomputedSamplePoints/pitch0.005_tolerance0',
              'precomputedSamplePoints/pitch0.005_tolerance0.005',
              'precomputedSamplePoints/pitch0.005_tolerance0.01',
              'precomputedSamplePoints/pitch0.005_tolerance0.02',
              'precomputedSamplePoints/pitch0.005_tolerance0.03',
              'precomputedSamplePoints/pitch0.005_tolerance0.04',
              'precomputedSamplePoints/pitch0.0075_tolerance0',
              'precomputedSamplePoints/pitch0.0075_tolerance0.005',
              'precomputedSamplePoints/pitch0.0075_tolerance0.01',
              'precomputedSamplePoints/pitch0.0075_tolerance0.02',
              'precomputedSamplePoints/pitch0.0075_tolerance0.03',
              'precomputedSamplePoints/pitch0.0075_tolerance0.04',
              'precomputedSamplePoints/pitch0.01_tolerance0',
              'precomputedSamplePoints/pitch0.01_tolerance0.005',
              'precomputedSamplePoints/pitch0.01_tolerance0.01',
              'precomputedSamplePoints/pitch0.01_tolerance0.02',
              'precomputedSamplePoints/pitch0.01_tolerance0.03',
              'precomputedSamplePoints/pitch0.01_tolerance0.04',
              'precomputedSamplePoints/pitch0.01_tolerance0.05',
              'precomputedSamplePoints/pitch0.02_tolerance0',
              'precomputedSamplePoints/pitch0.02_tolerance0.005',
              'precomputedSamplePoints/pitch0.02_tolerance0.01',
              'precomputedSamplePoints/pitch0.02_tolerance0.02',
              'precomputedSamplePoints/pitch0.02_tolerance0.03',
              'precomputedSamplePoints/pitch0.02_tolerance0.04',
              'precomputedSamplePoints/pitch0.02_tolerance0.05',
              'precomputedSamplePoints/pitch0.03_tolerance0',
              'precomputedSamplePoints/pitch0.03_tolerance0.005',
              'precomputedSamplePoints/pitch0.03_tolerance0.01',
              'precomputedSamplePoints/pitch0.03_tolerance0.02',
              'precomputedSamplePoints/pitch0.03_tolerance0.03',
              'precomputedSamplePoints/pitch0.03_tolerance0.04',
              'precomputedSamplePoints/pitch0.03_tolerance0.05',
              'precomputedSamplePoints/pitch0.04_tolerance0',
              'precomputedSamplePoints/pitch0.04_tolerance0.005',
              'precomputedSamplePoints/pitch0.04_tolerance0.01',
              'precomputedSamplePoints/pitch0.04_tolerance0.02',
              'precomputedSamplePoints/pitch0.04_tolerance0.03',
              'precomputedSamplePoints/pitch0.04_tolerance0.04',
              'precomputedSamplePoints/pitch0.04_tolerance0.05',
              'precomputedSamplePoints/pitch0.05_tolerance0',
              'precomputedSamplePoints/pitch0.05_tolerance0.005',
              'precomputedSamplePoints/pitch0.05_tolerance0.01',
              'precomputedSamplePoints/pitch0.05_tolerance0.02',
              'precomputedSamplePoints/pitch0.05_tolerance0.03',
              'precomputedSamplePoints/pitch0.05_tolerance0.04',
              'precomputedSamplePoints/pitch0.05_tolerance0.05',
              'precomputedSamplePoints/pitch0.06_tolerance0',
              'precomputedSamplePoints/pitch0.06_tolerance0.005',
              'precomputedSamplePoints/pitch0.06_tolerance0.01',
              'precomputedSamplePoints/pitch0.06_tolerance0.02',
              'precomputedSamplePoints/pitch0.06_tolerance0.03',
              'precomputedSamplePoints/pitch0.06_tolerance0.04',
              'precomputedSamplePoints/pitch0.06_tolerance0.05',
              'precomputedSamplePoints/pitch0.07_tolerance0',
              'precomputedSamplePoints/pitch0.07_tolerance0.005',
              'precomputedSamplePoints/pitch0.07_tolerance0.01',
              'precomputedSamplePoints/pitch0.07_tolerance0.02',
              'precomputedSamplePoints/pitch0.07_tolerance0.03',
              'precomputedSamplePoints/pitch0.07_tolerance0.04',
              'precomputedSamplePoints/pitch0.07_tolerance0.05',
              'precomputedSamplePoints/pitch0.08_tolerance0',
              'precomputedSamplePoints/pitch0.08_tolerance0.005',
              'precomputedSamplePoints/pitch0.08_tolerance0.01',
              'precomputedSamplePoints/pitch0.08_tolerance0.02',
              'precomputedSamplePoints/pitch0.08_tolerance0.03',
              'precomputedSamplePoints/pitch0.08_tolerance0.04',
              'precomputedSamplePoints/pitch0.08_tolerance0.05',
              'precomputedSamplePoints/pitch0.09_tolerance0',
              'precomputedSamplePoints/pitch0.09_tolerance0.005',
              'precomputedSamplePoints/pitch0.09_tolerance0.01',
              'precomputedSamplePoints/pitch0.09_tolerance0.02',
              'precomputedSamplePoints/pitch0.09_tolerance0.03',
              'precomputedSamplePoints/pitch0.09_tolerance0.04',
              'precomputedSamplePoints/pitch0.09_tolerance0.05',
              'precomputedSamplePoints/pitch0.12_tolerance0',
              'precomputedSamplePoints/pitch0.12_tolerance0.005',
              'precomputedSamplePoints/pitch0.12_tolerance0.01',
              'precomputedSamplePoints/pitch0.12_tolerance0.02',
              'precomputedSamplePoints/pitch0.12_tolerance0.03',
              'precomputedSamplePoints/pitch0.12_tolerance0.04',
              'precomputedSamplePoints/pitch0.12_tolerance0.05',
              'precomputedSamplePoints/pitch0.14_tolerance0',
              'precomputedSamplePoints/pitch0.14_tolerance0.005',
              'precomputedSamplePoints/pitch0.14_tolerance0.01',
              'precomputedSamplePoints/pitch0.14_tolerance0.02',
              'precomputedSamplePoints/pitch0.14_tolerance0.03',
              'precomputedSamplePoints/pitch0.14_tolerance0.04',
              'precomputedSamplePoints/pitch0.14_tolerance0.05',
              'precomputedSamplePoints/pitch0.16_tolerance0',
              'precomputedSamplePoints/pitch0.16_tolerance0.005',
              'precomputedSamplePoints/pitch0.16_tolerance0.01',
              'precomputedSamplePoints/pitch0.16_tolerance0.02',
              'precomputedSamplePoints/pitch0.16_tolerance0.03',
              'precomputedSamplePoints/pitch0.16_tolerance0.04',
              'precomputedSamplePoints/pitch0.16_tolerance0.05',
              'precomputedSamplePoints/pitch0.18_tolerance0',
              'precomputedSamplePoints/pitch0.18_tolerance0.005',
              'precomputedSamplePoints/pitch0.18_tolerance0.01',
              'precomputedSamplePoints/pitch0.18_tolerance0.02',
              'precomputedSamplePoints/pitch0.18_tolerance0.03',
              'precomputedSamplePoints/pitch0.18_tolerance0.04',
              'precomputedSamplePoints/pitch0.18_tolerance0.05',
              'precomputedSamplePoints/pitch0.1_tolerance0',
              'precomputedSamplePoints/pitch0.1_tolerance0.005',
              'precomputedSamplePoints/pitch0.1_tolerance0.01',
              'precomputedSamplePoints/pitch0.1_tolerance0.02',
              'precomputedSamplePoints/pitch0.1_tolerance0.03',
              'precomputedSamplePoints/pitch0.1_tolerance0.04',
              'precomputedSamplePoints/pitch0.1_tolerance0.05',
              'precomputedSamplePoints/pitch0.2_tolerance0',
              'precomputedSamplePoints/pitch0.2_tolerance0.005',
              'precomputedSamplePoints/pitch0.2_tolerance0.01',
              'precomputedSamplePoints/pitch0.2_tolerance0.02',
              'precomputedSamplePoints/pitch0.2_tolerance0.03',
              'precomputedSamplePoints/pitch0.2_tolerance0.04',
              'precomputedSamplePoints/pitch0.2_tolerance0.05',
              'precomputedSamplePoints/pitches',
              'precomputedSamplePoints/tolerances'],
              install_dir: precomputed_point_dir)

# pkg-config file
pkg = import('pkgconfig')
pkg.generate(xgandalf,
             filebase: 'xgandalf',
             extra_cflags: extra_cflag_list,
             description: 'Extended Gradient Descent Algorithm for Lattice Finding')

xgandalf_dep = declare_dependency(include_directories: incdir,
                                  link_with: xgandalf)

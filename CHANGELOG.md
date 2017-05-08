# UNANIMITY - CHANGELOG

## [3.1.0]

### Changed
 - Per ZMW timings are default on in DIAGNOSTICS mode or available via hidden
   option --zmwTimings. Output is BAM tag ms

## [3.0.0]

### Refactored
 - MultiMolecularIntegrator renamed to just Integrator
 - MonoMolecularIntegrator removed, all integrators now accept multiple molecules
 - VirtualTemplate removed, as without MonoMolecular it is no longer needed
 - MutatedTemplate added as a View object over some const template
 - Template::Mutate() now returns a MutatedTemplate instead of modifying the Template
 - Template was promoted from a member of Recursor to a member of EvaluatorImpl
 - Recursor refactored to take a template as an argument in most functions
 - Existing model files updated to match the new parent Recursor class
 - s/PB_CHEMISTRY_BUNDLE_DIR/SMRT_CHEMISTRY_BUNDLE_DIR/g

## [2.1.0]

### Added
 - Use pbcopper's q-gram index for sparse alignment
 - Replaced seqan MSA in ChimeraLabeler
 - support loading bundle models from PB_CHEMISTRY_BUNDLE_DIR
   environment variable

## [2.0.4]

### Added
 - Add pbcopper's ToolContract, summary is no longer a second output file
 - Differentiate between .xml and .bam output type
 - Enforce .pbi generation

## [2.0.3]

### Added
 - Switch from cpp-optparse to pbcopper, use pbcopper's CLI parsing

## [2.0.2]

### Added
 - Fix index errors in the Hirschberg aligner
 - Added a cleaner interface for AddRead/GetTemplate

## [2.0.1]

### Added
 - Add new ReleaseWithAssert CMAKE_BUILD_TYPE
 - Bump version (to cc2 + ccs)
 - Unify CCS and CC2 versioning under unanimity
 - Cleanup python/swig generation
 - Cleanup version handling

## [0.0.1]

### Added
 - Unify code base, refactor directory structure
 - Add pbccs, ConsensusCore2, pbsparse, and pbchimera
 - Code coverage report
 - Initial framework including pbbam, htslib, pbcopper

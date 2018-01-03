###########
# headers #
###########

if not meson.is_subproject()
  # pacbio/
  install_headers(
    files([
      'pacbio/UnanimityConfig.h',
      'pacbio/UnanimityVersion.h']),
    subdir : 'pacbio')

  # pacbio/align
  install_headers(
    files([
      'pacbio/align/AffineAlignment.h',
      'pacbio/align/AlignConfig.h',
      'pacbio/align/BandedAligner.h',
      'pacbio/align/BandedChainAlignment.h',
      'pacbio/align/ChainScorer.h',
      'pacbio/align/ChainSeeds.h',
      'pacbio/align/ChainSeedsConfig.h',
      'pacbio/align/FilterSeeds.h',
      'pacbio/align/FindSeeds.h',
      'pacbio/align/FindSeedsConfig.h',
      'pacbio/align/HomopolymerHasher.h',
      'pacbio/align/LinearAlignment.h',
      'pacbio/align/LocalAlignment.h',
      'pacbio/align/PairwiseAlignment.h',
      'pacbio/align/SeedScorer.h',
      'pacbio/align/SparseAlignment.h']),
    subdir : 'pacbio/align')

  # pacbio/align/internal
  install_headers(
    files([
      'pacbio/align/internal/BCAlignBlocks.h',
      'pacbio/align/internal/BCAlignImpl.h']),
    subdir : 'pacbio/align/internal')

  # pacbio/ccs
  install_headers(
    files([
      'pacbio/ccs/Consensus.h',
      'pacbio/ccs/ConsensusSettings.h',
      'pacbio/ccs/SparseAlignment.h',
      'pacbio/ccs/Whitelist.h',
      'pacbio/ccs/Consensus.h',
      'pacbio/ccs/ConsensusSettings.h',
      'pacbio/ccs/SparseAlignment.h',
      'pacbio/ccs/Whitelist.h']),
    subdir : 'pacbio/ccs')

  # pacbio/chimera
  install_headers(
    files([
      'pacbio/chimera/ChimeraLabel.h',
      'pacbio/chimera/ChimeraLabeler.h',
      'pacbio/chimera/ChimeraResultWriter.h']),
    subdir : 'pacbio/chimera')

  # pacbio/consensus
  install_headers(
    files([
      'pacbio/consensus/AbstractMatrix.h',
      'pacbio/consensus/Coverage.h',
      'pacbio/consensus/EasyReadScorer.h',
      'pacbio/consensus/Evaluator.h',
      'pacbio/consensus/Integrator.h',
      'pacbio/consensus/IntervalMask.h',
      'pacbio/consensus/MatrixViewConvention.h',
      'pacbio/consensus/ModelConfig.h',
      'pacbio/consensus/ModelSelection.h',
      'pacbio/consensus/Mutation.h',
      'pacbio/consensus/Polish.h',
      'pacbio/consensus/PolishResult.h',
      'pacbio/consensus/Template.h']),
    subdir : 'pacbio/consensus')

  # pacbio/consensus/internal
  install_headers(
    files([
      'pacbio/consensus/internal/ModelInternalInitializer.h']),
    subdir : 'pacbio/consensus/internal')

  # pacbio/data
  install_headers(
    files([
      'pacbio/data/ChemistryMapping.h',
      'pacbio/data/ChemistryTriple.h',
      'pacbio/data/Interval.h',
      'pacbio/data/IntervalTree.h',
      'pacbio/data/PlainOption.h',
      'pacbio/data/Read.h',
      'pacbio/data/ReadId.h',
      'pacbio/data/Sequence.h',
      'pacbio/data/State.h',
      'pacbio/data/StrandType.h',
      'pacbio/data/SubreadResultCounter.h']),
    subdir : 'pacbio/data')

  # pacbio/data/internal
  install_headers(
    files([
      'pacbio/data/internal/BaseEncoding.h',
      'pacbio/data/internal/ConversionFunctions.h']),
    subdir : 'pacbio/data/internal')

  # pacbio/denovo
  install_headers(
    files([
      'pacbio/denovo/PoaConsensus.h',
      'pacbio/denovo/PoaGraph.h',
      'pacbio/denovo/RangeFinder.h',
      'pacbio/denovo/SparsePoa.h']),
    subdir : 'pacbio/denovo')

  # pacbio/exception
  install_headers(
    files([
      'pacbio/exception/CCSExceptions.h',
      'pacbio/exception/InvalidEvaluatorException.h',
      'pacbio/exception/ModelError.h',
      'pacbio/exception/StateError.h']),
    subdir : 'pacbio/exception')

  # pacbio/io
  install_headers(
    files([
      'pacbio/io/Utility.h']),
    subdir : 'pacbio/io')

  # pacbio/parallel
  install_headers(
    files([
      'pacbio/parallel/WorkQueue.h']),
    subdir : 'pacbio/parallel')

  # pacbio/parallel
  install_headers(
    files([
      'pacbio/parallel/WorkQueue.h']),
    subdir : 'pacbio/parallel')

  # pacbio/util
  install_headers(
    files([
      'pacbio/util/ExecUtils.h',
      'pacbio/util/Timer.h']),
    subdir : 'pacbio/util')
endif

uny_include_directories = include_directories(['.', 'pacbio'])

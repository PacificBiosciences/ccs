###########
# headers #
###########

# pacbio/align
install_headers([
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
  'pacbio/align/SparseAlignment.h'],
  subdir : 'pacbio/align')

# pacbio/align/internal
install_headers([
  'pacbio/align/internal/BCAlignBlocks.h',
  'pacbio/align/internal/BCAlignImpl.h'],
  subdir : 'pacbio/align/internal')

# pacbio/ccs
install_headers([
  'pacbio/ccs/Consensus.h',
  'pacbio/ccs/ConsensusSettings.h',
  'pacbio/ccs/SparseAlignment.h',
  'pacbio/ccs/Whitelist.h',
  'pacbio/ccs/Consensus.h',
  'pacbio/ccs/ConsensusSettings.h',
  'pacbio/ccs/SparseAlignment.h',
  'pacbio/ccs/Whitelist.h'],
  subdir : 'pacbio/ccs')

# pacbio/chimera
install_headers([
  'pacbio/chimera/ChimeraLabel.h',
  'pacbio/chimera/ChimeraLabeler.h',
  'pacbio/chimera/ChimeraResultWriter.h'],
  subdir : 'pacbio/chimera')

# pacbio/consensus
install_headers([
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
  'pacbio/consensus/Template.h'],
  subdir : 'pacbio/consensus')

# pacbio/data
install_headers([
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
  'pacbio/data/SubreadResultCounter.h'],
  subdir : 'pacbio/data')

# pacbio/denovo
install_headers([
  'pacbio/denovo/PoaConsensus.h',
  'pacbio/denovo/PoaGraph.h',
  'pacbio/denovo/RangeFinder.h',
  'pacbio/denovo/SparsePoa.h'],
  subdir : 'pacbio/denovo')

# pacbio/exception
install_headers([
  'pacbio/exception/CCSExceptions.h',
  'pacbio/exception/InvalidEvaluatorException.h',
  'pacbio/exception/ModelError.h',
  'pacbio/exception/StateError.h'],
  subdir : 'pacbio/exception')

# pacbio/io
install_headers([
  'pacbio/io/Utility.h'],
  subdir : 'pacbio/io')

# pacbio/parallel
install_headers([
  'pacbio/parallel/WorkQueue.h'],
  subdir : 'pacbio/parallel')

# pacbio/parallel
install_headers([
  'pacbio/parallel/WorkQueue.h'],
  subdir : 'pacbio/parallel')

# pacbio/util
install_headers([
  'pacbio/util/ExecUtils.h',
  'pacbio/util/Timer.h'],
  subdir : 'pacbio/util')

uny_include_directories = include_directories(['.', 'pacbio'])

subdir('pacbio')

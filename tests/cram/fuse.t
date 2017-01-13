Make sure that consensus of aligned reads to identical reference stay identical
  $ $__PBTEST_FUSE_EXE $TESTDIR/../data/fuse/identical.bam -o $CRAMTMP/
  $ tail -n 1 $CRAMTMP/identical.cons | tr -d '\n' > $CRAMTMP/identical.seq
  $ tail -n 1 $TESTDIR/../data/fuse/identicalref.fasta> $CRAMTMP/identical.ref
  $ diff $CRAMTMP/identical.seq $CRAMTMP/identical.ref

Reference is missing one codon, restore it
  $ $__PBTEST_FUSE_EXE $TESTDIR/../data/fuse/singleCodonMissing.bam -o $CRAMTMP/
  $ tail -n 1 $CRAMTMP/singleCodonMissing.cons | tr -d '\n' > $CRAMTMP/singleCodonMissing.seq
  $ tail -n 1 $TESTDIR/../data/fuse/identicalref.fasta> $CRAMTMP/identical.ref
  $ diff $CRAMTMP/singleCodonMissing.seq $CRAMTMP/identical.ref

Reference is missing two codons too close together. Only restore the first
  $ $__PBTEST_FUSE_EXE $TESTDIR/../data/fuse/doubleCodonMissing.bam -o $CRAMTMP/
  $ tail -n 1 $CRAMTMP/doubleCodonMissing.cons | tr -d '\n' > $CRAMTMP/doubleCodonMissing.seq
  $ tail -n 1 $TESTDIR/../data/fuse/doubleCodonMissingTruth.fasta> $CRAMTMP/doubleCodonMissingTruth.ref
  $ diff $CRAMTMP/doubleCodonMissing.seq $CRAMTMP/doubleCodonMissingTruth.ref

Alignment has multiple gaps, remove them
  $ $__PBTEST_FUSE_EXE $TESTDIR/../data/fuse/gaps.bam -o $CRAMTMP/
  $ tail -n 1 $CRAMTMP/gaps.cons | tr -d '\n' > $CRAMTMP/gaps.seq
  $ tail -n 1 $TESTDIR/../data/fuse/identicalref.fasta> $CRAMTMP/identical.ref
  $ diff $CRAMTMP/gaps.seq $CRAMTMP/identical.ref
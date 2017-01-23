  $ $TESTDIR/../../../scripts/minorvariant/julietflow $TESTDIR/../../data/julietflow/flea99_97_1.bam $TESTDIR/../../data/julietflow/hxb2.fasta 2> $CRAMTMP/julietPerformance

Keep true positive rate of 1
  $ cut -f 1 -d' ' $CRAMTMP/julietPerformance
  1

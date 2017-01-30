Test mimosa 97-1-1-1 mix
  $ $TESTDIR/../../../scripts/minorvariant/julietflow $TESTDIR/../../data/julietflow/flea99_97_1.bam $TESTDIR/../../data/julietflow/hxb2.fasta 2> $CRAMTMP/julietPerformanceMimosa

Keep true positive rate of 1
  $ cut -f 1 -d' ' $CRAMTMP/julietPerformanceMimosa
  0.6

Test seabiscuit 96-1-1-1-1 mix
  $ $TESTDIR/../../../scripts/minorvariant/julietflow $TESTDIR/../../data/julietflow/sb99_96_1.bam $TESTDIR/../../data/julietflow/hxb2.fasta 2> $CRAMTMP/julietPerformanceSeabiscuit

Keep true positive rate of 1
  $ cut -f 1 -d' ' $CRAMTMP/julietPerformanceSeabiscuit
  1

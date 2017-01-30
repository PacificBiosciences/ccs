Make sure that the json output stays until we make adjustments to the algo or workflow
  $ $__PBTEST_JULIET_EXE $TESTDIR/../data/juliet_hiv_3000_96-1-1-1-1.align.bam -c "<HIV>" 2> $CRAMTMP/julietPerformance

Keep true positive rate of 1
  $ cut -f 1 -d' ' $CRAMTMP/julietPerformance
  1

Do not increase FP rate
  $ FP=$(cut -f 5 -d' ' $CRAMTMP/julietPerformance)
  $ echo $FP'<=1' | bc -l
  1

Test if we mistakenly overwrite

  $ touch exists.fq
  $ $__PBTEST_CCS_EXE $TESTDIR/../data/0passes.bam exists.fq
  OUTPUT: file already exists: 'exists.fq'
  [1]

  $ $__PBTEST_CCS_EXE $TESTDIR/../data/0passes.bam --force --minPasses=0 --minPredictedAccuracy=0.85 exists.fq
  $ if [ -s exists.fq ] ; then echo "dirty"; else echo " clean"; fi
  dirty

  $ rm exists.fq
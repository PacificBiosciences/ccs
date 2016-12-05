Make sure that the json output stays until we make adjustments to the algo or workflow
  $ $__PBTEST_JULIET_EXE $TESTDIR/../data/1000.bam -e FLEA_RQ99
  $ diff $TESTDIR/../data/1000.json 1000.json
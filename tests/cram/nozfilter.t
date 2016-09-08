Test to ensure the Z filtering can be disabled, the NaN parsing appears a bit finicky on some platforms

  $ $__PBTEST_CCS_EXE --zmws=111715 --minZScore=NaN $TESTDIR/../data/tiny.bam nofilter.fq
  /bin/sh: 2: --zmws=111715: not found
  [127]
  $ cat nofilter.fq
  cat: nofilter.fq: No such file or directory
  [1]

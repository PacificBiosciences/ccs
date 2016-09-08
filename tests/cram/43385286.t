Test a tiny collection of a few ZMWs, write to FASTQ for inspection

  $ $__PBTEST_CCS_EXE --minReadScore=0 --minSnr=0 --minZScore=nan $TESTDIR/../data/43385286.bam 43385286.fq
  /bin/sh: 2: --minReadScore=0: not found
  [127]
  $ cat 43385286.fq
  cat: 43385286.fq: No such file or directory
  [1]

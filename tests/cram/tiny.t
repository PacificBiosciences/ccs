Test a tiny collection of a few ZMWs, write to FASTQ for inspection

  $ $__PBTEST_CCS_EXE $TESTDIR/../data/tiny.bam tiny.fq
  /bin/sh: 2: /home/UNIXHOME/bbowman/jira/SAT-43/unanimity/tests/cram/../data/tiny.bam: Permission denied
  [126]
  $ cat tiny.fq
  cat: tiny.fq: No such file or directory
  [1]

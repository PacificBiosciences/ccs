
Test ccs on 100 zmws from the lexogen-SIRV dataset

  $ $__PBTEST_CCS_EXE --logLevel=DEBUG --minZScore -100 --maxDropFraction 0.8 $TESTDIR/../data/100zmws.bam test.fq
  >|> \d{8} \d{2}:\d{2}:\d{2}\.\d{3} -|- DEBUG      -|- main -|- [0-9,a-f,x]+|| -|- Found consensus models for: (P6-C4, S/P1-C1) (re)
  >|> \d{8} \d{2}:\d{2}:\d{2}\.\d{3} -|- DEBUG      -|- main -|- [0-9,a-f,x]+|| -|- Using consensus models for: (P6-C4) (re)

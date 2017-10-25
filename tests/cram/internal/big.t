This is a larger scale test run on a nightly basis.  It runs a good
chunk of a movie through CCS and pushes the results through ccscheck.
For now we just watch that the results don't change; in the future we
might like to replace that with some bulk checks on alignment identity,
for example.

Run CCS:

  $ DATADIR=/pbi/dept/consensus/testdata/unanimity-nightly

  $ ${__PBTEST_CCS_EXE} --logLevel=DEBUG --logFile=ccs.log --zmws 1-187412 $DATADIR/ds.subreadset.xml out.bam

Run ccscheck, check output:

  $ /pbi/dept/consensus/ccscheck/bin/ccscheck out.bam stats /pbi/dept/consensus/references/lambdaNEB.fasta
  $ sort -t, -n -k1,1 -k2,2 stats/zmws.csv > zmws.sorted.csv
  $ zdiff -NrU1 zmws.sorted.csv /pbi/dept/consensus/testdata/unanimity-nightly/big.sorted.to187412.csv.gz

Copy outputs to lastrun directory, for later inspection if need be

  $ cp out.bam* $DATADIR/lastrun/
  $ cp zmws.sorted.csv $DATADIR/lastrun/
  $ cp ccs.log ccs_report.txt $DATADIR/lastrun/

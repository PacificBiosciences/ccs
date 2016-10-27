#!/bin/bash

if [ "$(ps -p "$$" -o comm=)" != "bash" ]; then
    # Taken from http://unix-linux.questionfor.info/q_unix-linux-programming_85038.html
    bash "$0" "$@"
    exit "$?"
fi

echo "# DEPENDENCIES"
echo "## Load modules"
source /mnt/software/Modules/current/init/bash
module load gcc/5.3.0 python/2.7.9 zlib/1.2.5 graphviz samtools

source venv/bin/activate

echo "# TEST"
echo "## Running internal tests"
export __DATADIR="/pbi/dept/consensus/testdata/unanimity-nightly"
export __CCS_EXE="$(pwd -P)/build/ccs"
( mkdir -p tests/cram/internal )
( echo $'

run ccs

  $ mkdir results
  $ ${__CCS_EXE} --zmws 1-187412 --force ${__DATADIR}/ds.subreadset.xml ${__DATADIR}/lastrun/out.tmp.bam
  $ mv ${__DATADIR}/lastrun/out.tmp.bam ${__DATADIR}/lastrun/out.bam

run ccscheck, check output

  $ /pbi/dept/consensus/ccscheck/bin/ccscheck ${__DATADIR}/lastrun/out.bam stats /pbi/dept/consensus/references/lambdaNEB.fasta
  $ sort -t, -n -k1,1 -k2,2 stats/zmws.csv > ${__DATADIR}/lastrun/zmws.sorted.csv
  $ diff -NrU1 ${__DATADIR}/lastrun/zmws.sorted.csv /pbi/dept/consensus/testdata/unanimity-nightly/zmws.sorted.to187412.csv

' > tests/cram/internal/big.t )
( cram -v --xunit-file=uny-internal-cram.xml tests/cram/internal/*.t )

deactivate

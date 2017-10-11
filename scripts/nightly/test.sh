#!/bin/bash
set -euo pipefail

echo "#############################"
echo "# LOAD MODULES"
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load gcc/6.4.0
module load zlib/1.2.8
module load htslib/1.3.1

module load anaconda
PYTHONUSERBASE=$PWD/unyve
PATH=$PWD/unyve/bin:$PATH
export PATH PYTHONUSERBASE

echo "#############################"
echo "# INTERNAL TESTS"
export
__PBTEST_CCS_EXE="$(pwd -P)/build/ccs" scripts/cram --xunit-file=build/uny-internal-cram.xml --verbose tests/cram/internal/*.t

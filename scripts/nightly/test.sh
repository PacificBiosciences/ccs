#!/bin/bash
set -euo pipefail

echo "# DEPENDENCIES"
echo "## Load modules"
source /mnt/software/Modules/current/init/bash
module load gcc/4.9.2 python/2.7.9 zlib/1.2.5 graphviz samtools

set +u
source venv/bin/activate
set -u

echo "# TEST"
echo "## Running internal tests"
export
__PBTEST_CCS_EXE="$(pwd -P)/build/ccs" cram -v --xunit-file=uny-internal-cram.xml tests/cram/internal/*.t

set +u
deactivate
set -u

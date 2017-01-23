#!/bin/bash
set -euo pipefail

source /mnt/software/Modules/current/init/bash
module load cram/0.7 samtools

(cd artifacts/blasr/ && for i in *tgz; do tar --overwrite -x -f $i; done)
if [ ! -f `pwd`/artifacts/blasr/lib/libhdf5_cpp.so.12 ]
then
    ln -s /mnt/software/h/hdf5-tools/1.8.16/centos-7/lib/libhdf5_cpp.so `pwd`/artifacts/blasr/lib/libhdf5_cpp.so.12
fi

export PATH=`pwd`/artifacts/blasr/bin:`pwd`/artifacts/unanimity:$PATH
set +u
export LD_LIBRARY_PATH=/mnt/software/h/hdf5-tools/1.8.16/centos-7/lib:`pwd`/artifacts/blasr/lib:$LD_LIBRARY_PATH
set -u

wdir=$(pwd)
scripts/cram tests/cram/siv/julietflow.t --xunit-file=${wdir}/juliet-cram.xml
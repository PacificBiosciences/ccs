#!/bin/bash
set -euo pipefail

# Function definitions
GetBBRepo () {
    echo "## Clone $1"
    ( cd $3 && git clone ssh://git@bitbucket.nanofluidics.com:7999/$1/$2)
}

# Main script

echo "#############################"
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash
module load git gcc/4.9.2 python/2.7.9 cmake cram swig ccache virtualenv zlib/1.2.5 ninja boost
unset PKG_CONFIG_LIST

echo "#############################"
echo "# EXTERNAL DEPENDENCIES"
echo "## Create external dependency directory"
if [ -d _deps ] ; then rm -rf _deps ; fi
mkdir _deps
echo "## Create reverse external dependency directory"
if [ -d _rev_deps ] ; then rm -rf _rev_deps ; fi
mkdir _rev_deps

GetBBRepo sat GenomicConsensus _rev_deps
GetBBRepo sat ConsensusCore _deps
GetBBRepo sl pbcommand _deps
GetBBRepo sat pbcore _deps
GetBBRepo sat pacbiotestdata _deps

echo "## Fetch submodules"
git submodule update --init --remote

echo "#############################"
echo "# PRE-BUILD HOOK"
echo "## Check formatting"
./tools/check-formatting --all

echo "#############################"
echo "# VIRTUALENV"
echo "## Create missing virtualenv"
if [ ! -d unyve ] ; then /mnt/software/v/virtualenv/13.0.1/virtualenv.py unyve ; fi

echo "## Get into virtualenv"
set +u
source unyve/bin/activate
set -u

echo "## Install pip modules"
pip install --upgrade pip
pip install numpy cython h5py pysam cram nose jsonschema avro pytest
( cd _deps/pbcommand && pip install --no-deps . )
( cd _deps/pbcore && pip install --no-deps . )

echo "## Install pacbiotestdata"
( cd _deps/pacbiotestdata && git lfs pull && make python )

echo "#############################"
echo "# BUILD"
echo "## Create build directory "
if [ ! -d build ] ; then mkdir build ; fi

echo "## Build source"
( cd build &&\
  rm -rf * &&\
  cmake -DCMAKE_BUILD_TYPE=ReleaseWithAssert -GNinja .. )
( cd build && ninja )

echo "## pip install CC2"
CMAKE_BUILD_TYPE=ReleaseWithAssert CMAKE_COMMAND=cmake VERBOSE=1 pip install --verbose --upgrade --no-deps .

echo "## Install ConsensusCore"
( cd _deps/ConsensusCore && python setup.py install --boost=$BOOST_ROOT )

echo "## Install GC"
( cd _rev_deps/GenomicConsensus && pip install --upgrade --no-deps --verbose . )

echo "#############################"
echo "# TEST"
echo "## Unanimity tests"
( cd build && ninja check )

echo "## CC2 version test"
python -c "import ConsensusCore2 ; print ConsensusCore2.__version__"

echo "## Test CC2 via GC"
( cd _rev_deps/GenomicConsensus && make check )

set +u
deactivate
set -u

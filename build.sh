#!/bin/bash
set -euo pipefail

# Function definitions
GetBBRepo () {
    echo "## $1"
    if [ ! -d $2/$1 ]; then
        echo "### Clone"
        ( cd $2 && git clone ssh://git@bitbucket.nanofluidics.com:7999/sat/$1)
    else
        echo "### Update"
        ( cd $2/$1 && git pull)
    fi
}

GetGHRepo () {
    echo "## $1"
    if [ ! -d $2/$1 ]; then
        echo "### Clone"
        ( cd $2 && git clone https://github.com/PacificBiosciences/$1)
    else
        echo "### Update"
        ( cd $2/$1 && git pull)
    fi
}

# Main script

echo "#############################"
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash
module load git gcc/5.3.0 python/2.7.9 cmake cram swig ccache virtualenv zlib/1.2.5 ninja boost

echo "#############################"
echo "# EXTERNAL DEPENDENCIES"
echo "## Create external dependency directory"
if [ ! -d _deps ] ; then mkdir _deps ; fi
echo "## Create reverse external dependency directory"
if [ ! -d _rev_deps ] ; then mkdir _rev_deps ; fi

GetBBRepo GenomicConsensus _rev_deps
GetBBRepo ConsensusCore _deps
GetGHRepo pbcommand _deps
GetGHRepo pbcore _deps
GetGHRepo PacBioTestData _deps

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
pip install numpy cython h5py pysam cram nose jsonschema avro
( cd _deps/pbcommand && pip install --no-deps . )
( cd _deps/pbcore && pip install --no-deps . )

echo "## Install PacBioTestData"
( cd _deps/PacBioTestData && git lfs pull && make python )

echo "#############################"
echo "# BUILD"
echo "## Create build directory "
if [ ! -d build ] ; then mkdir build ; fi

echo "## Build source"
( cd build &&\
  rm -rf * &&\
  CMAKE_BUILD_TYPE=ReleaseWithAssert cmake -DZLIB_INCLUDE_DIR=/mnt/software/z/zlib/1.2.5/include -DZLIB_LIBRARY=/mnt/software/z/zlib/1.2.5/lib/libz.so -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++" -GNinja .. )
( cd build && ninja htslibSrc )
( cd build && ninja )

echo "## pip install CC2"
CMAKE_BUILD_TYPE=ReleaseWithAssert CMAKE_COMMAND=cmake ZLIB_INCLUDE_DIR=/mnt/software/z/zlib/1.2.5/include ZLIB_LIBRARY=/mnt/software/z/zlib/1.2.5/lib/libz.so VERBOSE=1 pip install --verbose --upgrade --no-deps .

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
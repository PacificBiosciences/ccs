#!/bin/bash

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

echo "## GenomicConsensus"
if [ ! -d _rev_deps/GenomicConsensus ]; then
    echo "### Clone"
    ( cd _rev_deps && git clone ssh://git@bitbucket.nanofluidics.com:7999/sat/GenomicConsensus)
else
    echo "### Update"
    ( cd _rev_deps/GenomicConsensus && git pull)
fi

echo "## ConsensusCore"
if [ ! -d _deps/ConsensusCore ]; then
    echo "### Clone"
    ( cd _deps && git clone ssh://git@bitbucket.nanofluidics.com:7999/sat/ConsensusCore)
else
    echo "### Update"
    ( cd _deps/ConsensusCore && git pull)
fi

echo "## Fetch submodules"
git submodule update --init --remote --depth 1

echo "#############################"
echo "# PRE-BUILD HOOK"
echo "## Check formatting"
./tools/check-formatting --all

echo "#############################"
echo "# VIRTUALENV"
echo "## Create missing virtualenv"
if [ ! -d unyve ]
then 
    /mnt/software/v/virtualenv/13.0.1/virtualenv.py unyve
fi

echo "## Get into virtualenv"
source unyve/bin/activate

echo "## Install pip modules"
pip install --upgrade pip
pip install numpy cython h5py pysam cram nose jsonschema avro
pip install --no-deps git+https://github.com/PacificBiosciences/pbcommand.git
pip install --no-deps git+https://github.com/PacificBiosciences/pbcore.git

echo "## Get PacBioTestData"
if [ -d _deps/PacBioTestData ]
then
    rm -rf _deps/PacBioTestData
fi
( cd _deps                                                         &&\
git clone https://github.com/PacificBiosciences/PacBioTestData.git &&\
cd PacBioTestData                                                  &&\
git lfs pull                                                       &&\
make python )

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

deactivate
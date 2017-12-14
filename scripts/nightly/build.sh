#!/bin/bash
set -euo pipefail

echo "# DEPENDENCIES"
echo "## Load modules"
source /mnt/software/Modules/current/init/bash
module load git gcc python/2 cmake swig ccache zlib ninja boost

echo "## Get into virtualenv"
if [ ! -d venv ]
then
    virtualenv venv
fi
set +u
source venv/bin/activate
set -u

echo "## Install pip modules"
pip install --upgrade pip
pip install cram nose

pip install xmlbuilder jsonschema avro requests iso8601
pip install --no-deps git+https://github.com/PacificBiosciences/pbcommand.git

pip install numpy cython h5py pysam
pip install --no-deps git+https://github.com/PacificBiosciences/pbcore.git

# pip install jinja2 networkx xmlbuilder requests fabric
# pip install --no-deps git+https://github.com/PacificBiosciences/pbsmrtpipe.git

#echo "## Fetch unanimity submodules"
#( git submodule update --init --remote )

echo "# BUILD"
echo "## Create build directory "
if [ ! -d build ] ; then mkdir build ; fi

echo "## Build source"
( cd build &&\
  rm -rf * &&\
  CMAKE_BUILD_TYPE=ReleaseWithAssert cmake -DJULIET_INHOUSE_PERFORMANCE=T -DZLIB_INCLUDE_DIR=/mnt/software/z/zlib/1.2.5/include -DZLIB_LIBRARY=/mnt/software/z/zlib/1.2.5/lib/libz.so -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++" -GNinja .. )
( cd build && ninja )

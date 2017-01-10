#!/bin/bash
set -euo pipefail

echo "# DEPENDENCIES"
echo "## Load modules"
source /mnt/software/Modules/current/init/bash
module load git gcc/5.3.0 python/2.7.9 cmake swig ccache virtualenv zlib/1.2.5 ninja boost

echo "## Get into virtualenv"
if [ ! -d venv ]
then
    /mnt/software/v/virtualenv/13.0.1/virtualenv.py venv
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

echo "## Fetch unanimity submodules"
( git submodule update --init --remote )

echo "# BUILD"
echo "## Create build directory "
if [ ! -d build ] ; then mkdir build ; fi

echo "## Build source"
( cd build &&\
  rm -rf * &&\
  CMAKE_BUILD_TYPE=ReleaseWithAssert cmake -DJULIET_INHOUSE_PERFORMANCE=T -DZLIB_INCLUDE_DIR=/mnt/software/z/zlib/1.2.5/include -DZLIB_LIBRARY=/mnt/software/z/zlib/1.2.5/lib/libz.so -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++" -GNinja .. )
( cd build && ninja )

#!/bin/bash
set -euo pipefail

echo "#############################"
echo "# LOAD MODULES"
source /mnt/software/Modules/current/init/bash
module load git gcc/4.9.2 python/2.7.9 cmake cram swig ccache virtualenv zlib/1.2.5 ninja boost lcov gcovr

echo "## Fetch submodules"
git submodule update --init --remote

echo "#############################"
echo "# PRE-BUILD HOOK"
echo "## Check formatting"
./tools/check-formatting --all

echo "#############################"
echo "# BUILD"
echo "## Create build directory "
if [ ! -d build ] ; then mkdir build ; fi

echo "## Build source"
( cd build &&\
  rm -rf * &&\
  cmake -DCMAKE_BUILD_TYPE=Debug -DUNY_inc_coverage=ON -GNinja .. )
( cd build && ninja ccs test_unanimity )

echo "#############################"
echo "# TEST"
echo "## Unanimity tests"
( cd build && ./tests/test_unanimity && \
  ./ccs --zmws 109700 ../tests/data/tiny.bam /dev/null/tiny.fq && \
  find . -type f -iname '*.o' | xargs gcov -acbfu {} \; >/dev/null && \
  mkdir coverage && pushd coverage && mv ../*.gcov . && \
  sed -i -e 's@Source:@Source:../@' *.gcov && \
  sed -i -e 's@Graph:@Graph:../@' *.gcov && \
  sed -i -e 's@Data:@Data:../@' *.gcov )

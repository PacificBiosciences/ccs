#!/bin/bash
set -euo pipefail

# Function definitions
GetBBRepo () {
    type module >& /dev/null || . /mnt/software/Modules/current/init/bash
    module load git
    echo "## Clone $1"
    ( cd $3 && git clone ssh://git@bitbucket.nanofluidics.com:7999/$1/$2)
}

# Main script

echo "#############################"
echo "# LOAD MODULES"
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
module load gcc
module load ccache
module load boost
module load zlib
module load htslib
module load cmake
module load swig
module load ninja
if [[ $USER == "bamboo" ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.mobs.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi
export CCACHE_COMPILERCHECK='%compiler% -dumpversion'
export CCACHE_BASEDIR=$PWD

echo "#############################"
echo "# EXTERNAL DEPENDENCIES"
echo "## Create external dependency directory"
if [[ -d _deps ]] ; then rm -rf _deps ; fi
mkdir _deps
echo "## Create reverse external dependency directory"
if [[ -d _rev_deps ]] ; then rm -rf _rev_deps ; fi
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
rm -rf unyve
mkdir unyve
module load anaconda
PYTHONUSERBASE=$PWD/unyve
PATH=$PWD/unyve/bin:$PATH
export PATH PYTHONUSERBASE
NEXUS_WHEEL=http://nexus/repository/unsupported/pitchfork/gcc-6.4.0
PIP="pip --cache-dir=$PWD/.pip --disable-pip-version-check"

echo "## Install pip modules"
$PIP install --user \
  cram==0.7
$PIP install --user --no-index \
  $NEXUS_WHEEL/pythonpkgs/pysam-0.13-cp27-cp27mu-linux_x86_64.whl \
  $NEXUS_WHEEL/pythonpkgs/xmlbuilder-1.0-cp27-none-any.whl \
  $NEXUS_WHEEL/pythonpkgs/avro-1.7.7-cp27-none-any.whl \
  $NEXUS_WHEEL/pythonpkgs/iso8601-0.1.12-py2.py3-none-any.whl \
  $NEXUS_WHEEL/pythonpkgs/tabulate-0.7.5-cp27-none-any.whl \
  http://nexus/repository/unsupported/distfiles/coverage-4.4.1.tar.gz
$PIP install --user --no-index -e _deps/pbcommand
$PIP install --user --no-index -e _deps/pbcore

echo "## Install pacbiotestdata"
ln -sfn ../data _deps/pacbiotestdata/pbtestdata/data
$PIP install --user -e _deps/pacbiotestdata

echo "#############################"
echo "# BUILD"
echo "## Create build directory "
if [[ ! -d build ]] ; then mkdir build ; fi

echo "## Build source"
( cd build &&\
  rm -rf * &&\
  cmake -DCMAKE_BUILD_TYPE=ReleaseWithAssert -GNinja .. )
( cd build && sed -i -e 's@/-I/mnt/software@/ -I/mnt/software@g' build.ninja && ninja )

echo "## pip install CC2"
LDFLAGS='-static-libstdc++' \
CMAKE_BUILD_TYPE=ReleaseWithAssert \
CMAKE_COMMAND=cmake \
VERBOSE=1 \
$PIP install --user --no-index --verbose .
ldd unyve/lib/python2.7/site-packages/_ConsensusCore2.so

echo "## CC2 version test"
python -c "import ConsensusCore2 ; print ConsensusCore2.__version__"

rm -f unyve/bin/g++
module unload ccache
cat > g++ <<EOF
#!/bin/bash
$(which g++) -static-libstdc++ \$@
EOF
chmod +x g++
module load ccache
mv g++ unyve/bin/
echo "## Install ConsensusCore"
( cd _deps/ConsensusCore \
  && python setup.py \
       bdist_wheel \
       --boost=$BOOST_ROOT \
  && echo dist/ConsensusCore-*.whl | \
     xargs $PIP install --user --verbose )
rm unyve/bin/g++
ldd unyve/lib/python2.7/site-packages/_ConsensusCore.so

echo "## Install GC"
$PIP install --user --no-index --verbose _rev_deps/GenomicConsensus

echo "#############################"
echo "# TEST"
echo "## Unanimity tests"
( cd build && ninja check )

#echo "## Test CC2 via GC"
#module load cram
#module add mummer/3.23
#module add exonerate/2.0.0
#( cd _rev_deps/GenomicConsensus && make check )

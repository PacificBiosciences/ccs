
Building/installing pbccs
=========================

First make sure you have boost (>=1.5.8) and cmake (3.3) installed.

Next, here's a recipe for setting up the dependencies in the
_deps folder:

  ```sh
  $ mkdir _deps
  $ pushd _deps
  $ git clone https://github.com/PacificBiosciences/htslib.git
  $ git clone https://github.com/PacificBiosciences/pbbam.git
  $ git clone https://github.com/PacificBiosciences/ConsensusCore2.git
  $ git clone https://github.com/PacificBiosciences/seqan.git
  $ popd
  ```

Next, configure and build as below.  Sorry about the awkward call to
readlink below, we are working on fixing that.  Note that Mac users
should use GNU readlink (available as greadlink if you `brew install
coreutils`).

  ```sh
  $ mkdir build; pushd build
  $ cmake                                               \
      -DSEQAN_INCLUDE_DIRS=../_deps/seqan/include       \
      -DPacBioConsensus_RootDir=../_deps/ConsensusCore2 \
      -DPacBioBAM_RootDir=../_deps/pbbam                \
      -DHTSLIB_ROOTDIR=`readlink -f ../_deps/htslib`    \
      ..
  $ make
  ```


Building/installing bax2bam
===========================

If you want to use pbccs with data from the PacBio RSII instrument,
you will need the bax2bam converter.

(TODO: instructions on how to build it!!)

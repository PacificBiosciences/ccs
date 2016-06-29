*Note: all are welcome to build our software from source, but official
 support is only provided for official builds provided by PacBio
 (the SMRT Analysis suite)*

## Building/installing pbccs

### Easiest way: use `pitchfork`

The easiest way to build is to use the `pitchfork` build tool, which
automates dependency fetching/resolution:

  ```sh
  $ git clone https://github.com/PacificBiosciences/pitchfork
  $ cd pitchfork
  $ make pbccs
  ```

### Easy way: manually

Building from scratch requires system-wide installed boost (>=1.5.8), 
cmake (3.3), and a c++11 compiler (>=gcc-4.8, clang):

  ```sh
  $ git clone https://github.com/PacificBiosciences/pbccs && cd pbccs
  $ git submodule update --init --remote
  $ mkdir build && cd build
  $ cmake .. && make -j
  $ cd ..
  $ bin/ccs --help
  ```

## Building/installing bax2bam

If you want to use pbccs with data from the PacBio RSII instrument,
you will need the bax2bam converter, which you can build using
`pitchfork`.  At present, you need to build BLASR, as bax2bam is part
of the BLASR suite.

  ```sh
  $ git clone https://github.com/PacificBiosciences/pitchfork
  $ cd pitchfork
  $ make blasr
  ```

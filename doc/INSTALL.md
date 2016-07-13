*Note: all are welcome to build our software from source, but official
 support is only provided for official builds provided by PacBio
 (the SMRT Analysis suite)*

 ***

### Individual tools: use `pitchfork`

The easiest way to build is to use the `pitchfork` build tool, which
automates dependency fetching/resolution:

  ```sh
  git clone https://github.com/PacificBiosciences/pitchfork && \
  cd pitchfork                                              && \
  make pbccs
  ```

### All tools: manually

Building from scratch requires system-wide installed boost (>=1.5.8), 
cmake (3.2), and a c++11 compiler (>=gcc-5.3.0, clang):

  ```sh
  git clone https://github.com/PacificBiosciences/unanimity && \
  cd unanimity                                              && \
  git submodule update --init --remote                      && \
  mkdir build                                               && \
  cd build                                                  && \
  cmake -DCMAKE_INSTALL_PREFIX=~/bin ..                     && \
  make -j install                                           && \
  ~/bin/ccs
  ```

### Individual tools: manually

Invoke the different `make` targets, currently available

  ```sh
  $ make pbccs
  ```
*Note: all are welcome to build our software from source, but official
 support is only provided for official builds provided by PacBio
 (the SMRT Analysis suite)*

 ***

# There isn't a proper way for installation.

### Install ConsensusCore2 python library for GenomicConsensus

Building from scratch requires system-wide installed boost (>=1.5.8), 
cmake (3.2), python 2.x, and a c++11 compiler (>=gcc-5.3.0, clang):

  ```sh
  git clone https://github.com/PacificBiosciences/unanimity                               && \
  cd unanimity                                                                            && \
  git submodule update --init --remote                                                    && \
  pip install numpy cython h5py pysam cram nose                                           && \
  pip install --upgrade --no-deps git+https://github.com/PacificBiosciences/pbcommand.git && \
  pip install --upgrade --no-deps git+https://github.com/PacificBiosciences/pbcore.git    && \
  pip install --upgrade --no-deps .
  ```

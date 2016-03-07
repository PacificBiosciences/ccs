
## Installation on the PacBio cluster:

First, set up the build environment.

```sh
module purge   # start with a clean environment
module load gcc/4.8.4 boost/1.58 swig/2.0.10 cmake/3.2.2
```

Second: are you building `ConsensusCore2` for use by
`GenomicConsensus` (Python; requires `ConsensusCore2` SWIG bindings)
or `pbccs`/`pblaa` (pure C++)?


### Build just the C++ library

```sh
mkdir build; cd build; cmake ..
make
```

(In fact you probably don't have to do this directly; building
`pbccs`/`pblaa` should do this as a subtask)


### Build the Python SWIG bindings

(Note: this builds the C++ library automatically; you don't need to
build it as a separate step)

We recommend installing ConsensusCore2 within a Python `virtualenv`.

    ```sh
    python setup.py install
    ```

<p align="center">
  <img src="doc/img/unanimity.png" alt="unanimity logo"/>
</p>
<h1 align="center">Unanimity</h1>
<p align="center">C++ library and its applications to generate and process accurate consensus sequences</p>

***
## Documentation

 - [Getting Started](doc/INSTALL.md)
 - Projects
   - Available
     - Consensus Core
     - [Circular Consensus Calling `ccs`](doc/PBCCS.md)
   - Work in Progress
     - Genomic Consensus Calling `gcpp`
     - Viral Haplotype Phasing `eden`
   - Planned
     - Minor Variant Calling
 - [Developer environment](doc/DEVELOPER.md)
 - [PacBio open source license](LICENSE)

## Quick Tools Overview

### [Circular Consensus Calling](doc/PBCCS.md)

`ccs` takes multiple reads of the same SMRTbell sequence and combines 
them, employing a statistical model, to produce one high quality consensus sequence.
    
### [Genomic Consensus Calling](doc/GCPP.md)

`gcpp` will replace the current python [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus), until then please use the existing solution.

### Minor Variant Calling

This tool will be able to precisely call single-nucleotide variants from consensus data.

### [Viral Haplotype Phasing](doc/EDEN.md)

`eden` will leverage CCS reads to identify low-frequency haplotypes within polyploid samples.

## Help

Support is only provided for official and stable
[SMRT Analysis builds](http://www.pacb.com/products-and-services/analytical-software/)
provided by PacBio and not for source builds.
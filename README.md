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
     - [Minor Variant Calling `juliet`](doc/JULIET.md)
   - Planned
     - Genomic Consensus Calling `gcpp`
     - Viral Haplotype Phasing `eden`
 - [Developer environment](doc/DEVELOPER.md)
 - [PacBio open source license](LICENSE)

## Quick Tools Overview

### [Circular Consensus Calling](doc/PBCCS.md)

`ccs` takes multiple reads of the same SMRTbell sequence and combines
them, employing a statistical model, to produce one high quality consensus sequence.

### Genomic Consensus Calling

`gcpp` will replace the current python [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus), until then please use the existing solution.

### [Minor variant caller](doc/JULIET.md)

`juliet` identifies minor variants from aligned ccs reads.

### Viral Haplotype Phasing

`eden` will leverage CCS reads to identify low-frequency haplotypes within polyploid samples.

## Help

Support is only provided for official and stable
[SMRT Analysis builds](http://www.pacb.com/products-and-services/analytical-software/)
provided by PacBio and not for source builds.
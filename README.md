<p align="center">
  <img src="doc/img/unanimity.png" alt="unanimity logo"/>
</p>
<h1 align="center">Unanimity</h1>
<p align="center">C++ library and its applications to generate and process accurate consensus sequences</p>

***

## Tools

### Circular Consensus Calling

`ccs` takes multiple reads of the same SMRTbell sequence and combines 
them, employing a statistical model, to produce one high quality consensus sequence.
    
### Genomic Consensus Calling

`gcpp` will replace the current python [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus), until then please use the existing solution.

### Minor Variant Calling

This tool will be able to precisely call single-nucleotide variants from consensus data.

### Viral Haplotype Phasing

`eden` will leverage CCS reads to identify low-frequency haplotypes within polyploid samples.

## Content

 - [Installation](doc/INSTALL.md)
 - Tools and Algorithms
   - [pbccs](doc/PBCCS.md)
   - [gcpp](doc/GCPP.md)
   - [eden](doc/EDEN.md)

## Help

Issues? Bugs? Please create a github issue.

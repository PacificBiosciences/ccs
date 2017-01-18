<h1 align="center">
    fuse - Reduce alignment into its representative sequence
</h1>

<p align="center">
  <img src="img/fuse.png" alt="Logo of Fuse" width="100px"/>
</p>

## Install
Install the unanimity suite and one of the binaries is called `fuse`.

## Input data
*Fuse* operates on aligned records in the BAM format.
BAM files has to PacBio-compliant, meaning, cigar `M` is forbidden.

## Scope
Current scope of *Fuse* is creation of a high-quality consensus sequence.
Fuse includes in-frame insertions with a certain distance to each other.
Major deletions are being removed.

## Output
*Fuse* provides a FASTA file per input. If no output prefix has been specified
via `-o`, the output files are called `inputPrefix.cons`, otherwise
`outputPrefix_inputPrefix.cons`. There is no warning, existing files will be
overwritten.

## Example
Simple example:
```
fuse m530526.align.bam
```

Output: `m530526.cons`
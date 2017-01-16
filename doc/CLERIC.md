<h1 align="center">
    Cleric - Swap BAM alignment reference
</h1>

<p align="center">
  <img src="img/cleric.png" alt="Logo of Cleric" width="200px"/>
</p>

## Install
Install the unanimity suite and one of the binaries is called `cleric`.

## Input data
*Cleric* operates on aligned records in the BAM format, the original reference
and the target reference as FASTA.
BAM file has to PacBio-compliant, meaning, cigar `M` is forbidden.
Two sequences have to be provided, either in individual files or combined in one.
The header of the original reference must match the reference name in the BAM.

## Scope
Current scope of *Cleric* is converting a given alignment to a different
reference. This is done by aligning the original and target reference sequences.
A transitive alignment is used to generate the new alignment.

## Output
*Cleric* provides a BAM file with the file named as provided via `-o`.

## Example
Simple example:
```
cleric m530526.align.bam reference.fasta new_ref.fasta
```

Or:
```
cat reference.fasta new_ref.fasta > combined.fasta
cleric m530526.align.bam combined.fasta
```
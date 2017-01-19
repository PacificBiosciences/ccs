<h1 align="center">
    julietflow - Minor variant pipeline
</h1>

## Install
The script is located under `scripts/minovariant/julietflow`.

## Input data
*Julietflow* operates on unaligned ccs reads in the BAM format and a close
reference sequence.
BAM files has to PacBio-compliant, meaning, cigar `M` is forbidden.

## Scope
Current scope of *julietflow* is to automatize the re-align workflow for
targeted HIV POL minor variant calling.

## Output
*Julietflow* provides the output of juliet in the current directory and the
intermediate files in tmp/.

## Example
```
julietflow m530526.ccs.bam hxb2.fasta
```

Output: `m530526_cleric.html`

## Dependencies
*Julietflow* currently depends on `juliet`, `cleric`, `fuse` all of which are
part of unanimity, `blasr`, and `samtools`. Those are all part of the
SMRT-bundle that is officially being provided by PacBio.

## Workflow
<p align="center">
  <img src="img/julietflow.png" alt="Julietflow workflow" width="500px"/>
</p>
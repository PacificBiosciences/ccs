<h1 align="center">
    juliet - Minor Variant Caller
</h1>

<p align="center">
  <img src="img/juliet.png" alt="Logo of Juliet" width="400px"/>
</p>

## Install
Install the unanimity suite and one of the binaries is called `juliet`.

## Input data
*Juliet* operates on aligned CCS records in the BAM format.
Reads should be created with [CCS2](../PBCCS.md).
The reference has to be [HXB2](https://www.hiv.lanl.gov/components/sequence/HIV/asearch/query_one.comp?se_id=K03455).
BAM files have to PacBio-compliant, meaning, cigar `M` is forbidden.
*Juliet* does not demultiplex barcoded data, provide one BAM per barcode.

## Scope
Current scope of *Juliet* is identification of single-nucleotide variants in the
function region POL of HIV.

## Output
*Juliet* identifies SNVs and provides for each gene
(*Protease*, *Reverse Transcriptase*, *RNase*, *Integrase*)
mutated amino acids. For each amino acid, each nucleotide of the codon is
being reported along with the observed relative abundance and its p-value;
for major substitutions the p-value is replaced by `M`. Known
positions of drug resistance mutations are marked.

Example:

```
Reverse Transcriptase
=====================
  41 K(AAA) => N[(A 1 M)](A 0.0801 2.02769e-24)](C 0.3594 3.89532e-97)] <+> NNRTI surveillance
```

At RT amino acid position 41, the HXB2 reference is Lysine `K` with codon `AAA`.
The observed amino acid is Asparagine `N` with codon `AAC`.
The first nucleotide is wild-type `A`,
the second `A` is only present as a minor variant with an observed frequency
of 8.01% and a p-value of 2e-24,
and the third `C` is a high-frequent variant of 35.94% and a p-value of 3e-97.

## Usage

Typical usage is:
```
juliet -i path/to/aligned_ccs_input.bam
```

## Additional output
*Juliet* produces a `msa_output.tsv` file, containing positional counts and
their p-values.
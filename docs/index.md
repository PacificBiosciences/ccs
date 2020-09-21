---
layout: default
title: CCS Home
nav_order: 1
description: "Generate Highly Accurate Single-Molecule Consensus Reads (HiFi Reads)."
permalink: /
---

<p align="center">
  <img src="img/ccs_card.png" alt="CCS logo" width="650px"/>
</p>

***

_ccs_ combines multiple subreads of the same SMRTbell molecule using a
statistical model to produce one highly accurate consensus sequence,
also called a HiFi read, along with base quality values.
This tool powers the _Circular Consensus Sequencing_ workflow in SMRT Link.

## Availability
The latest `ccs` can be installed via the bioconda package `pbccs`.

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
for information on Installation, Support, License, Copyright, and Disclaimer.

## Latest Version
Version **5.0.0**: [Full changelog here](/changelog)

## Schematic Workflow
<p align="center"><img width="1000px" src="img/generate-hifi.png"/></p>

## Execution
**Input**: Subreads from a single movie in PacBio BAM format (`.subreads.bam`).

**Output**: Consensus reads in a format inferred from the file extension:
unaligned BAM (`.bam`); bgzipped FASTQ (`.fastq.gz`);
or SMRT Link XML (`.consensusreadset.xml`) which also generates a corresponding
unaligned BAM file.

Run on a full movie:

    ccs movie.subreads.bam movie.ccs.bam

Parallelize by using `--chunk`.
See [how-to chunk](/faq/parallelize).

Feel free to [make _ccs_ slightly faster via environment variables](/faq/performance#can-i-tune-performance-without-sacrificing-output-quality).

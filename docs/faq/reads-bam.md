---
layout: default
parent: FAQ
title: reads.bam
---

# What is the `reads.bam`?
Have you ever run _ccs_ with different cutoffs, e.g. tuning `--min-rq` , because
out of the fear of missing out on yield?
Similar to the CLR instrument mode, in which subreads are accompanied by
a scraps file, _ccs_ offers a new mode to never lose a single read due to
filtering, without massive run time increase by polishing low-pass productive ZMWs.

Starting with SMRT Link v10.0 and Sequel IIe, _ccs_ v5.0 or newer is able to generate
one representative sequence per productive ZMW, irrespective of quality and passes.
This ensures no yield loss due to filtering and enables users to have maximum
control over their data. Never fear again that SMRT Link or the Sequel IIe
HiFi mode filtered precious data.

**Attention:** If you work with the `reads.bam` file directly, be aware that CCS reads of all
qualities are present. This file needs to be understood before piping
into your typical HiFi application.

## How to generate `reads.bam`?

The default command-line behavior has not changed;
it still generates only HiFi quality reads by default.
But the new `--all` mode has been set as default when running the
_Circular Consensensus Sequencing_ SMRT Link application or
selecting the on-instrument Sequel IIe capabilities:
<p align="left"><img width="500px" src="../img/run-design-oiccs.png"/></p>

## What is in the `reads.bam`?

- HiFi Reads with predicted accuracy ≥Q20 (`rq ≥ 0.99`)
- Lower-quality but still polished consensus reads with predicted accuracy <Q20 (`rq < 0.99`)
- Unpolished consensus reads (`rq = -1`)
- Partial or single full-length subreads unaltered (`rq = -1`)

## How to get HiFi reads

### SMRT Link

If you want to only use HiFi reads, SMRT Link automatically generates additional
files for your convenience that only contain HiFi reads:

 - hifi_reads.**fastq**.gz
 - hifi_reads.**fasta**.gz
 - hifi_reads.**bam**

### Command line

Following tools can be installed with

    conda install -c bioconda tool_name

#### extracthifi
We provide a simple tool, called `extracthifi` to generate a HiFi-only BAM from a `reads.bam` file. Usage is:

    extracthifi reads.bam extracthifi.bam

#### bamtools
Alternatively use `bamtools`:

    bamtools filter -in reads.bam -out hifi_reads.bam -tag "rq":">=0.99"

## FAQ: How can I filter by number of passes?

We **strongly** advise against filtering by anything than predicted accuracy,
BAM tag `rq`. The `rq` tag is the best predictor for read quality. Number of
passes is not reliable enough and you might discard too much data. This `np`
tag is an implementation detail that is not guaranteed to be present in future
_ccs_ versions.

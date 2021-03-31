---
layout: default
parent: FAQ
title: SQIIe
---

# Sequel IIe (SQIIe) FAQ
## Is _ccs_ the same for on-instrument Sequel IIe and off-instrument SMRT Link?
Yes. The source code is identical for a given release, e.g. ICS v10.0 and SMRT
Link v10.0. The only difference is how we generate the _ccs_ binaries for
*SQIIe* and *SMRT Link*; each binary is compiled for different environments,
with varying degrees of hardware specializations and optimizations. _ccs_
binaries for SQIIe are slightly faster than for SMRT Link, but all binaries from
a given release yield identical results.

## Will it be like this in the future?
There are no plans to bifurcate. In fact, we dedicate a lot of energy into
having all x86 binaries compile from the same source code and produce identical
results.

## What is the difference between bioconda and official PacBio releases?
Bioconda releases are only available as developer versions. They may contain
features and bug fixes that have not been integrated in our release software yet.


## Which _ccs_ options are used for on-instrument processing?
SQIIe default: `--all`

SQIIe with kinetic information: `--all --all-kinetics --subread-fallback`

## But the `reads.bam` header contains significantly more options!
Yes, those are mostly used to explicitly name auxilliary files. Everything explained:

```
/opt/pacbio/pa-ccs/current/bin/ccs \ # binary
    /data/pa/m64012_201202_224049.consensusreadset.xml \ # XML output name
    --all --all-kinetics --subread-fallback\ # all mode with kinetics
    --streamed \ # do not use a BAM file as input, but a BAM input stream
    --log-level INFO \
    --suppress-reports \ # do not generate report or metric files per default
    --log-file /data/pa/m64012_201202_224049.ccs.log \ # log file name
    --bam /data/pa/m64012_201202_224049.reads.bam \ # explicit BAM file name
    --report-json /data/pa/m64012_201202_224049.ccs_reports.json \ # the ccs_reports file as JSON
    --report-file /data/pa/m64012_201202_224049.ccs_reports.txt \ # the ccs_reports file as human readable
    --metrics-json /data/pa/m64012_201202_224049.zmw_metrics.json.gz \ # metrics file name
    --hifi-summary-json /data/pa/m64012_201202_224049.hifi_summary.json \ # summary JSON file for hifi statistics
    --stderr-json-log # instrument-specific additional JSON logging
```



## What files (formats) are transferred from the instruments?
For HiFi runs, the main file is a `reads.bam` file, prefixed by the movie name, and it is accompanied by two
core [auxilliary files](/faq/reports-aux-files):
* `*.ccs_reports.log`
* `*.zmw_metrics.json.gz`

## What is in the `reads.bam`?
The on-instrument _ccs_ version and also SMRT Link ≥v10 run in the `--all` mode
by default. In this mode, _ccs_ outputs one representative sequence per
productive ZMW, irrespective of quality and passes. More information
[in the `--all` mode FAQ](/faq/mode-all) and [in the `reads.bam` FAQ](/faq/reads-bam).

## Can you go back to `subreads.bam` from `reads.bam`?
Not when operating the instrument in CCS mode. See next question.

## Can I get subreads and CCS reads transferred for the same SMRT Cell?
You can select either CLR or CCS reads coming off the instrument. If you want
CLR and CCS reads, you can select HiFi generation on SMRT Link, which
transfers CLR off instrument and runs CCS on your SMRT Link instance:
<p align="left"><img width="800px" src="../img/run-design-slccs.png"/></p>

## What is the expected run time for _ccs_ on the instrument and will it delay the next SMRT cell?
The design of SQIIe allows to have _ccs_ running in parallel with the
acquisition of the subsequent SMRT cell, and will not interfere or delay subsequent acquisitions.
The hardware and software is configured to accommodate future increases in HiFi yield.

## How does _ccs_ on- and off-instrument performance compare?
Performance of _ccs_ on the instrument is roughly equivalent to running on a
dedicated system with 384 logical cores. On-instrument processing includes
parallel post-primary computation and result transfers up to 90% faster than for
CLR runs.

## Can you specify custom parameters on the instrument, like by-strand CCS?
There is no possibility to specify custom parameters. We are open to feedback
what might be of general interest in future instrument software updates.

## Can you include kinetics information or is it included by design?
You may select to retain kinetic information in Run-Design:
<p align="left"><img width="500px" src="../img/run-design-kinetics.png"/></p>
Be aware that the output BAM files will more than 5x larger due to additional
kinetics.

## Can I generate a `reads.bam` from a `subreads.bam` file which was generated with older chemistry / instrument software (even from the Sequel System)?
Generating a `reads.bam` has been possible since _ccs_ v5 and supported
chemistries can be found [here](/faq/chemistry).

## What is the file size of a HiFi dataset transferred from the instruments for one SMRT Cell on the SQIIe System?
We provide upper bounds and average file sizes how the BAM file size scales with
HiFi yield in the [BAM output FAQ](/faq/bam-output.html#how-does-the-output-bam-file-size-scale-with-yield).
From this page, a typical 30 GBases HiFi yield run will result in a ~33 GByte
`reads.bam` file with default parameters and ~174 GB with HiFi kinetics
included.

## Can demultiplexing of barcoded datasets be done on the SQIIe instrument?
No, this currently not supported. Demultiplexing continues to be supported through SMRT Link or with lima on the command line.

## Can the `reads.bam` still be used in a CLR application? E.g. to make full use of all the data, not just the HiFi reads?
In SMRT Link and in general, only the HiFi reads are used from the `reads.bam`.
It should not be necessary to include other reads. Some pipelines do allow
including CCS reads <Q20 (e.g. Iso-Seq). You are free to use the all reads from
the `reads.bam` files with your own tools.

## Can I extract the different subtypes of reads from the new `reads.bam`?
Yes, the SMRT Analysis pipeline "Export Reads" in SMRT Link v10.0 or newer can
export HiFi reads to BAM/FASTA/FASTQ format; when adjusting minimum CCS
predicted accuracy, you can include CCS reads <Q20. On the command line, tools
can be used to filter the BAM file for the read quality `rq` tag, please see
the [`reads.bam` FAQ](/faq/reads-bam).

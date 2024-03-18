---
layout: default
parent: FAQ
title: Fail reads
---

# What is the fail_reads.bam file?
On the Revio platform, for each `hifi_reads.bam` there is one `fail_reads.bam`
file. This file contains one representative read per ZMW. The `ff` tag encodes
why each read ended up in this file:

|  Flag  |                                                       Description                                                        |
| ------ | ------------------------------------------------------------------------------------------------------------------------ |
| `0x1`  | CCS reads with predicted accuracy below QV 20.                                                                           |
| `0x2`  | Control CCS reads.                                                                                                       |
| `0x4`  | Single-stranded CCS reads.                                                                                               |
| `0x8`  | The median full-length subread from molecules that do not produce a CCS read but have at least one full pass.            |
| `0x10` | CCS reads which are a concatenation of the adapter, with possible short non-adapter sequence in between.                 |
| `0x20` | CCS reads with miscalled adapter which is enclosed by a sequence and its reverse complement, either spanning to the end. |
| `0x40` | CCS reads that have one or more adapters close to either end.                                                            |


---
layout: default
parent: FAQ
title: BAM output
---

## What BAM tags are generated?

|  Tag  | Type  |                                                                                                    Description                                                                                                    |
| :---: | :---: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ac`  | `B,i` | [Detected and missing adapter counts](/faq/missing-adapters)                                                                                                                                                      |
| `ec`  |  `f`  | [Effective coverage](/faq/accuracy-vs-passes#how-is-number-of-passes-computed)                                                                                                                                    |
| `fi`  | `B,C` | [Double-strand forward IPD (codec V1)](/faq/kinetics)                                                                                                                                                             |
| `ff`  |  `i`  | [Fail reads](/faq/fail-reads)                                                                                                                                                                                     |
| `fn`  |  `i`  | [Double-strand forward number of complete passes (zero or more)](/faq/kinetics)                                                                                                                                   |
| `fp`  | `B,C` | [Double-strand forward PulseWidth (codec V1)](/faq/kinetics)                                                                                                                                                      |
| `ip`  | `B,C` | [Single-strand IPD (codec V1)](/faq/kinetics)                                                                                                                                                                     |
| `ma`  |  `i`  | [Missing adapters bitmask](/faq/missing-adapters)                                                                                                                                                                 |
| `np`  |  `i`  | [Number of full-length subreads](/faq/accuracy-vs-passes#how-is-number-of-passes-computed)                                                                                                                        |
| `pw`  | `B,C` | [Single-strand PulseWidth (codec V1)](/faq/kinetics)                                                                                                                                                              |
| `ri`  | `B,C` | [Double-strand reverse IPD (codec V1)](/faq/kinetics)                                                                                                                                                             |
| `rn`  |  `i`  | [Double-strand reverse number of complete passes (zero or more)](/faq/kinetics)                                                                                                                                   |
| `rp`  | `B,C` | [Double-strand reverse PulseWidth (codec V1)](/faq/kinetics)                                                                                                                                                      |
| `rq`  |  `f`  | [Predicted average read accuracy](/how-does-ccs-work#9-qv-calculation)                                                                                                                                            |
| `sa`  | `B,I` | [Run-length encoded per-base coverage by subread alignments in form of <length>,<coverage>,...](https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html#use-of-read-tags-for-hifi-per-read-base-pileup-summary) |
| `sm`  | `B,C` | [Per-base number of aligned matches](https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html#use-of-read-tags-for-hifi-per-read-base-pileup-summary)                                                            |
| `sx`  | `B,C` | [Per-base number of aligned mismatches](https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html#use-of-read-tags-for-hifi-per-read-base-pileup-summary)                                                         |
| `sn`  | `B,f` | Signal-to-noise ratios for each nucleotide                                                                                                                                                                        |
| `zm`  |  `i`  | ZMW hole number                                                                                                                                                                                                   |
| `RG`  |  `z`  | Read group                                                                                                                                                                                                        |


## How does the output BAM file size scale with yield?
For each base, the output BAM file size scales as follows
 - 0.5 byte/base for the actual base (4-bit encoding)
 - 1 byte/base for the QV
 - 1 byte/base for the forward PW
 - 1 byte/base for the forward IPD
 - 1 byte/base for the reverse PW
 - 1 byte/base for the reverse IPD

For a normal _ccs_ run without kinetics, the upper bound is 1.5 bytes/base.
If _ccs_ is run **with** kinetics, the upper bound is 5.5 bytes/base.

Per-read meta information add a fixed amount of 32 bytes per read:
 - `ec`,`rq` : float, each 4 bytes
 - `sn`: float array, 4x4 bytes
 - `np`, `zm`: int32_t, 4 byte
 - `RG`: string of length 8, 8x1 bytes

The actual output BAM that _ccs_ generates is compressed. Compression is
data-dependent and because of that, upper bounds can't be provided.

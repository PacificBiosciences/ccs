---
layout: default
parent: FAQ
title: BAM output
---

## What BAM tags are generated?

|  Tag  | Type  | Description |
| :---: | :---: | ----------- |
| `ec`  | `f`   | [Effective coverage](/faq/accuracy-vs-passes#how-is-number-of-passes-computed)|
| `fi`  | `B,C` | [Forward IPD (codec V1)](/faq/kinetics)|
| `fn`  | `i`   | [Forward number of complete passes (zero or more)](/faq/kinetics)|
| `fp`  | `B,C` | [Forward PulseWidth (codec V1)](/faq/kinetics)|
| `np`  | `i`   | [Number of full-length subreads](/faq/accuracy-vs-passes#how-is-number-of-passes-computed)|
| `ri`  | `B,C` | [Reverse IPD (codec V1)](/faq/kinetics)|
| `rn`  | `i`   | [Reverse number of complete passes (zero or more)](/faq/kinetics)|
| `rp`  | `B,C` | [Reverse PulseWidth (codec V1)](/faq/kinetics)|
| `rq`  | `f`   | [Predicted average read accuracy](/how-does-ccs-work#9-qv-calculation)|
| `sn`  | `B,f` | Signal-to-noise ratios for each nucleotide|
| `zm`  | `i`   | ZMW hole number |
| `RG`  | `z`   | Read group |


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
For a 19kb insert library and 30h movie time, the _ccs_ BAM files scale on
average with:

| BAM name             |        Options                             | Bytes/<br>Base | Bytes/<br>HiFiBase | Example<br>(GBytes) | Example<br>(GBytes) |
| -------------------- | ------------------------------------------ | :------------: | :----------------: | :-----------------: | :-----------------: |
| hifi.bam             |                                            | 0.7            | 0.7                | 100                 | 63                  |
| hifi.hifikin.bam     | `--hifi-kinetics`                          | 3.7            | 3.7                | 528                 | 336                 |
| reads.bam            | `--all`                                    | 0.55           | 1.1                | 157                 | 100                 |
| reads.hifikin.bam    | `--all --hifi-kinetics`                    | 2.3            | 4.5                | 642                 | 409                 |
| reads.allkin.bam     | `--all --all-kinetics`                     | 2.9            | 5.7                | 814                 | 518                 |
| reads.allkin.sub.bam | `--all --all-kinetics --subread-fallback`  | 3.0            | 5.8                | 828                 | 527                 |

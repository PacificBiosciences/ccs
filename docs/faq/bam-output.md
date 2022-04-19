---
layout: default
parent: FAQ
title: BAM output
---

## What BAM tags are generated?

|  Tag  | Type  |                                        Description                                         |
| :---: | :---: | ------------------------------------------------------------------------------------------ |
| `ac`  | `B,i` | [Detected and missing adapter counts](/faq/missing-adapters)                               |
| `ec`  |  `f`  | [Effective coverage](/faq/accuracy-vs-passes#how-is-number-of-passes-computed)             |
| `fi`  | `B,C` | [Double-strand forward IPD (codec V1)](/faq/kinetics)                                      |
| `fn`  |  `i`  | [Double-strand forward number of complete passes (zero or more)](/faq/kinetics)            |
| `fp`  | `B,C` | [Double-strand forward PulseWidth (codec V1)](/faq/kinetics)                               |
| `ip`  | `B,C` | [Single-strand IPD (codec V1)](/faq/kinetics)                                              |
| `ma`  |  `i`  | [Missing adapters bitmask](/faq/missing-adapters)                                          |
| `np`  |  `i`  | [Number of full-length subreads](/faq/accuracy-vs-passes#how-is-number-of-passes-computed) |
| `pw`  | `B,C` | [Single-strand PulseWidth (codec V1)](/faq/kinetics)                                       |
| `ri`  | `B,C` | [Double-strand reverse IPD (codec V1)](/faq/kinetics)                                      |
| `rn`  |  `i`  | [Double-strand reverse number of complete passes (zero or more)](/faq/kinetics)            |
| `rp`  | `B,C` | [Double-strand reverse PulseWidth (codec V1)](/faq/kinetics)                               |
| `rq`  |  `f`  | [Predicted average read accuracy](/how-does-ccs-work#9-qv-calculation)                     |
| `sn`  | `B,f` | Signal-to-noise ratios for each nucleotide                                                 |
| `zm`  |  `i`  | ZMW hole number                                                                            |
| `RG`  |  `z`  | Read group                                                                                 |


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

| Read types                 | Kinetics                 | Options                                   | Bytes/<br>Base | Bytes/<br>HiFiBase | Example<br>(GBytes) | Example<br>(GBytes) |
|----------------------------|--------------------------|-------------------------------------------|:--------------:|:------------------:|:-------------------:|:-------------------:|
| HiFi                       | None                     |                                           |      0.7       |        0.7         |         100         |         63          |
| HiFi                       | HiFi                     | `--hifi-kinetics`                         |      3.7       |        3.7         |         528         |         336         |
| HiFi + LQ CCS + unpolished | None                     | `--all`                                   |      0.55      |        1.1         |         157         |         100         |
| HiFi + LQ CCS + unpolished | HiFi                     | `--all --hifi-kinetics`                   |      2.3       |        4.5         |         642         |         409         |
| HiFi + LQ CCS + unpolished | HiFi + LQ CCS            | `--all --all-kinetics`                    |      2.9       |        5.7         |         814         |         518         |
| HiFi + LQ CCS + fallback   | HiFi + LQ CCS + fallback | `--all --all-kinetics --subread-fallback` |      3.0       |        5.8         |         828         |         527         |

**Legend:**
 - `HiFi` - Polished CCS reads with predicted accuracy greater equals Q20, optionally with kinetics
 - `LQ CCS` - Polished CCS reads with predicted accuracy below Q20, optionally with kinetics
 - `unpolished` - Unpolished consensus sequence with two or fewer passes, no kinetics possible
 - `fallback` - One representative subread for ZMWs, instead of an unpolished consensus sequence, optionally with kinetics

The Sequel IIe system either runs with `--all` per default or optionally with `--all --all-kinetics --subread-fallback`.

---
layout: default
parent: FAQ
title: Mode --split-heteroduplex
---

# Attention: This is an early access feature!

## What is heteroduplex filtering?
Starting with _ccs_ v6.2.0, single-strand artifacts, such as insertions larger
than 20 bases, do not necessarily have to be filtered out. Using `--split-heteroduplexes`,
_ccs_ is able to split ZMW on-the-fly after detecting such heteroduplex and
process each strand separately. As a consequence, _ccs_ has to distinguish between
double-stranded (DS) and single-stranded (DS) ZMWs and their consensus reads.
Implications:

 * Single-strand reads are stored in an extra file
 * Summary logs report double-strand and single-strand metrics
 * `ccs_reports.txt` file contains two columns, double-strand and single-strand reads

We are currently investigating how reliable we can detect indel and SNV
heteroduplexes and might add those to strand-aware splitting in future versions.

### Additional `*.stranded.bam` file
The file `outputPrefix.stranded.bam` contains all single-strand reads. Read names
follow the by-strand scheme with `/fwd` and `/rev` suffixes. There are up to two
reads per split ZMW.

### Summary logs
At the end of each execution, _ccs_ reports for `--log-level INFO` a summary.
This summary contains combined and individual metrics for DS and SS.

```
-------------------------------------------------
Summary stats abbreviations:
ZMW         - A productive Zero-Mode Waveguide
DS          - Double Strand
SS          - Single Strand
DS-ZMW      - All subreads were used from a single ZMW
SS-ZMW      - ZMW is split into fwd and rev strands,
              each strand is polished individually
DS-Read     - CCS read of a DS-ZMW
SS-Read     - CCS read of one strand of a SS-ZMW
HiFi        - CCS reads with predicted accuracy >=Q20
UMY         - Unique Molecular Yield of all reads passing filters
HiFi Yield  - UMY of >=Q20 DS- and SS-ZMWs, longest read per ZMW
-------------------------------------------------
ZMWs Input    : 53895
ZMWs Written  : 22684
 - DS / SS    : 22644 / 40
UMY           : 413.2 MBases (6.8 GBases/hr)
 - DS / SS    : 412.4 MBases / 733.7 KBases
HiFi Yield    : 413.5 MBases (6.8 GBases/hr)
 - DS / SS    : 412.4 MBases / 1.0 MBases
HiFi Reads    : 22701
 - DS / SS    : 22644 / 57
HiFi Avg Size : 18.2 KBases
HiFi Avg QV   : 30.2
```

### Strand-aware `ccs_reports.txt`
Typical content of the strand-aware `ccs_reports.txt` file. Contrary to the
default output, this file does not report numbers in ZMWs, but actual DS and SS
reads. Accounting in SS ZMWs is not possible, as one strand might fail and the
other succeed.

```
                           Double-Strand Reads  Single-Strand Reads
Inputs                   :      53590 (99.43%)         609 (0.564%)

Passed                   :      22644 (42.25%)          57 (9.360%)
Failed                   :      30946 (57.75%)         552 (90.64%)

Tandem repeats           :        461 (1.490%)           0 (0.000%)

Exclusive failed counts
Below SNR threshold      :        870 (2.811%)           0 (0.000%)
Median length filter     :          0 (0.000%)           0 (0.000%)
Shortcut filters         :          0 (0.000%)           0 (0.000%)
Lacking full passes      :      26226 (84.75%)           0 (0.000%)
Coverage drops           :         30 (0.097%)           0 (0.000%)
Insufficient draft cov   :         61 (0.197%)         310 (56.16%)
Draft too different      :          0 (0.000%)           0 (0.000%)
Draft generation error   :        173 (0.559%)          54 (9.783%)
Draft above --max-length :          0 (0.000%)           0 (0.000%)
Draft below --min-length :          0 (0.000%)           0 (0.000%)
Reads failed polishing   :          0 (0.000%)           0 (0.000%)
Empty coverage windows   :          3 (0.010%)           0 (0.000%)
CCS did not converge     :          2 (0.006%)           0 (0.000%)
CCS below minimum RQ     :       3581 (11.57%)         188 (34.06%)
Unknown error            :          0 (0.000%)           0 (0.000%)
```

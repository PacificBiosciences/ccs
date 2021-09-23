---
layout: default
parent: FAQ
title: Mode --by-strand
---

## Can I produce one consensus sequence for each strand of a molecule?
Yes, please use `--by-strand`. Make sure that you have sufficient coverage,
as `--min-passes` are per strand in this case. For each strand, _ccs_
generates one consensus read that has to pass all filters.
Read name suffix indicates strand. Example:

    m64011_190714_120746/14/ccs/rev
    m64011_190714_120746/35/ccs/fwd

How does `--by-strand` work? For each ZMW:
 * Determine orientation of reads with respect to the one closest to the median length
 * Sort reads into two buckets, forward and reverse strands
 * Treat each strand as an individual entity as we do with ZMWs
   * Apply all filters per strand individually
   * Create a draft for each strand
   * Polish each strand
   * Write out each polished strand consensus

### Summary logs
At the end of each execution, _ccs_ reports a summary for `--log-level INFO`.
This summary contains combined and individual metrics for double- and single-strand reads.
Using `--by-strand`, only single-strand reads are reported.

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
ZMWs Input    : 5390
ZMWs Written  : 1061
 - DS / SS    : 0 / 1061
UMY           : 18.8 MBases (2.2 GBases/hr)
 - DS / SS    : 0 Bases / 18.8 MBases
HiFi Yield    : 33.7 MBases (3.9 GBases/hr)
 - DS / SS    : 0 Bases / 33.7 MBases
HiFi Reads    : 1909
 - DS / SS    : 0 / 1909
HiFi Avg Size : 17.7 KBases
HiFi Avg QV   : 23.0
```

### By-strand `ccs_reports.txt`
Typical content of the by-strand `ccs_reports.txt` file. Contrary to the
default output, this file does not report numbers in ZMWs, but single-strand
reads. Accounting in single-strand ZMWs is not possible, as one strand might fail
and the other succeed.

```
                           Single-Strand Reads
Inputs                   :       9313 (86.38%)

Passed                   :       1909 (20.50%)
Failed                   :       7404 (79.50%)

Tandem repeats           :         69 (0.932%)

Exclusive failed counts
Below SNR threshold      :        101 (1.364%)
Median length filter     :          0 (0.000%)
Shortcut filters         :          0 (0.000%)
Lacking full passes      :       4888 (66.02%)
Coverage drops           :          6 (0.081%)
Insufficient draft cov   :         30 (0.405%)
Draft too different      :          0 (0.000%)
Draft generation error   :         39 (0.527%)
Draft above --max-length :          0 (0.000%)
Draft below --min-length :          0 (0.000%)
Reads failed polishing   :          0 (0.000%)
Empty coverage windows   :          1 (0.014%)
CCS did not converge     :          0 (0.000%)
CCS below minimum RQ     :       2339 (31.59%)
Unknown error            :          0 (0.000%)
```

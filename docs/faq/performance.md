---
layout: default
parent: FAQ
title: Performance
---

## How fast is CCS?
We tested CCS runtime using 500 ZMWs per length bin with exactly 7 passes.

<img width="1000px" src="../img/runtime.png"/>

### How does that translate into time to result per SMRT Cell?
We will measure time to result for Sequel System and Sequel System II CCS sequencing collections
on a PacBio recommended HPC, according to the
[Sequel II System Compute Requirements](https://www.pacb.com/wp-content/uploads/SMRT_Link_Installation_v701.pdf)
with 192 physical or 384 hyper-threaded cores.

1) Sequel System: 15 kb insert size, 24-hours movie, 37 GB raw yield, 2.3 GB HiFi UMY
2) Sequel II System: 15 kb insert size, 30-hours movie, 340 GB raw yield, 24 GB HiFi UMY

| CCS version | Sequel System | Sequel II System |
| :-: | :-: | :-: |
| ≤3.0.0 | 1 day | >1 week |
| 3.4.1 | 3 hours | >1 day |
| 4.0.0 | 40 minutes | 6 hours |
| ≥4.2.0 | **30 minutes** | **4 hours** |

### How is CCS speed affected by raw base yield?
Raw base yield is the sum of all polymerase read lengths.
A polymerase read consists of all subreads concatenated
with SMRTbell adapters in between.

Raw base yield can be increased with
1) higher percentage of single-loaded ZMWs and
2) longer movie times that lead to longer polymerase read lengths.

Since the first version, _ccs_ scales linear in (1) the number of single loaded
ZMWs per SMRT Cell.
Starting with version 3.3.0 _ccs_ scales linear in (2) the polymerase read length
and with version 4.0.0 _ccs_ scales sublinear.

### What did change in each version?

| CCS version | O(insert size) | O(#passes) |
| :-: | :-: | :-: |
| ≤3.0.0 | quadratic | linear |
| 3.4.1 | **linear** | linear |
| ≥4.0.0 | linear | **sublinear** |

### How can version 4.0.0 be sublinear in the number of passes?
With the introduction of improved heuristics, individual draft bases can skip
polishing if they are of sufficient quality.
The more passes a ZMW has, the fewer bases need additional polishing.

## Can I tune _ccs_ to get improved results?
No, we optimized _ccs_ such that there is a good balance between speed and
output quality.

## Does speed impact quality and yield?
Yes it does. With ~35x speed improvements from version 3.1.0 to 4.0.0 and
consequently reducing CPU time from >60,000 to <2,000 core hours,
heuristics and changes in algorithms lead to slightly lower yield and
accuracy if run head-to-head on the same data set. Internal tests show
that _ccs_ 4.0.0 introduces no regressions in CCS-only Structural Variant
calling and has minimal impact on SNV and indel calling in DeepVariant.
In contrast, lower DNA quality has a bigger impact on quality and yield.

## Can I tune performance without sacrificing output quality?
The bioconda _ccs_ ≥v5.0 binaries statically link [mimalloc](https://github.com/microsoft/mimalloc).
Depending on your system, additional performance tuning can be achieved.
Internally, we use following mimalloc environment variables to improve _ccs_
performance.

```
MIMALLOC_PAGE_RESET=0 MIMALLOC_LARGE_OS_PAGES=1 ccs <movie>.subreads.bam ccs.bam --log-level INFO
```

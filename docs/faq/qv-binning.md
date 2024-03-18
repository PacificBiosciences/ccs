---
layout: default
parent: FAQ
title: File compression
---

## File compression

File compression is important to increase storage density of genomes per hard
drive and reduce time to transfer files. We offer two approaches to
significantly reduce file size.

### QV binning

Compression of per-base quality values (QVs) is an effective method to reduce
file size. We adopt per-base QV binning in _ccs_ and can achieve up to 40%
reduction in BAM file size with similar SNV and InDel detection performance.

**How it works**: QV compression is performed on base QVs after read quality
`rq` has been computed. For each consensus read, per-base QVs are assigned to
seven bins and fixed average QVs are assigned to each bin

|    QV bin    | Mean QV | ASCII |
| ------------ | ------- | ----- |
| `[  0,  6 ]` | Q3      | `$`   |
| `[  7, 13 ]` | Q10     | `+`   |
| `[ 14, 19 ]` | Q17     | `2`   |
| `[ 20, 24 ]` | Q22     | `7`   |
| `[ 25, 29 ]` | Q27     | `<`   |
| `[ 30, 39 ]` | Q35     | `D`   |
| `[ 40, 93 ]` | Q40     | `I`   |

### Sorting by sequence similarity

In addition to QV filtering, sorting unaligned HiFi BAM files by sequence similarity
can reduce file size by up to an additional 30% using

    samtools sort -M -K 20 hifi_reads.bam -o hifi_reads.sorted.bam



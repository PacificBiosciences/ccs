---
layout: default
parent: FAQ
title: Parallelize
---

## How can I parallelize on multiple servers?
You can parallelize by chunking. Since _ccs_ v4.0.0, direct chunking via `--chunk`
is possible. For this, the `.subreads.bam` file must accompanied by a
`.pbi` file. To generate the index `subreads.bam.pbi`, use
`pbindex`, which can be installed with `conda install pbbam`.

    pbindex movie.subreads.bam

An example workflow, all ccs invocations can run simultaneously:

    ccs movie.subreads.bam movie.ccs.1.bam --chunk 1/10 -j <THREADS>
    ccs movie.subreads.bam movie.ccs.2.bam --chunk 2/10 -j <THREADS>
    ...
    ccs movie.subreads.bam movie.ccs.10.bam --chunk 10/10 -j <THREADS>

Merge chunks with `pbmerge` and index with `pbindex`

    pbmerge -o movie.ccs.bam movie.ccs.*.bam
    pbindex movie.ccs.bam

or use `samtools`

    samtools merge -@8 movie.ccs.bam movie.ccs.*.bam

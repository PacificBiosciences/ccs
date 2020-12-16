---
layout: default
parent: FAQ
title: Low complexity
---

## Does _ccs_ dislike low-complexity regions?
Low-complexity comes in many shapes and forms.
A particular challenge for _ccs_ are highly enriched tandem repeats, like
hundreds of copies of `AGGGGT`.
Prior _ccs_ v5.0, inserts with many copies of a small repeat likely not generate
a consensus sequence.
Since _ccs_ v5.0, every ZMW is tested if it contains a tandem repeat
of length `--min-tandem-repeat-length 1000`.
For this, we use [symmetric DUST](https://doi.org/10.1089/cmb.2006.13.1028)
and in particular the [sdust](https://github.com/lh3/sdust) implementation,
but slightly modified.
If a ZMW is flagged as a tandem repeat, internally `--disable-heuristics`
is activated for only this ZMW, and various filters that are known to exclude
low-complexity sequences are disabled.
This recovers most of the low-complexity consensus sequences, without impacting
run time performance.

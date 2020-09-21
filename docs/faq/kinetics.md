---
layout: default
parent: FAQ
title: Kinetics
---

## Is it possible to use HiFi reads to call base modifications?
Base modifications can be inferred from per-base pulse width (PW) and
inter-pulse duration (IPD) kinetics.
Running _ccs_ with `--hifi-kinetics` generates averaged kinetic information
for polished reads, independently for both strands of the insert.
Forward is defined with respect to the orientation represented in ``SEQ`` and
is considered to be the native orientation. As with other PacBio-specific
tags, aligners will not re-orient these fields.

Minor cases exist where a certain orientation may get filtered out entirely
from a ZMW, preventing valid values from being passed for that record. In
these cases, empty lists will be passed for the respective record/orientation
and number of passes will be set to zero.

In order to facilitate the use of HiFi reads with base modifications workflows,
we have added an executable in pbbam called `ccs-kinetics-bystrandify` which
creates a pseudo `--by-strand` BAM with corresponding `pw` and `ip` tags
that imitates a normal, unaligned subreads BAM. You can install pbbam from
Bioconda by calling `conda install pbbam`.

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

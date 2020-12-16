---
layout: default
parent: FAQ
title: Mode --all
---

# Process _all_ reads
Have you ever run _ccs_ with different cutoffs, e.g. tuning `--min-rq` , because
out of the fear of missing out on yield?
Similar to the CLR instrument mode, in which subreads are accompanied by
a scraps file, _ccs_ offers a new mode to never lose a single read due to
filtering, without massive run time increase by polishing low-pass productive ZMWs.

Starting with SMRT Link v10.0 and Sequel IIe, _ccs_ v5.0 or newer is able to generate
one representative sequence per productive ZMW, irrespective of quality and passes.
This ensures no yield loss due to filtering and enables users to have maximum
control over their data. Never fear again that SMRT Link or the Sequel IIe
HiFi mode filtered precious data.

The default command-line behavior has not changed;
it still generates only HiFi quality reads by default.
But the new `--all` mode has been set as default when running the
_Circular Consensensus Sequencing_ SMRT Link application or
selecting the on-instrument Sequel IIe capabilities:
<p align="left"><img width="500px" src="../img/run-design-oiccs.png"/></p>

Output will be a `reads.bam` that contains:

- HiFi Reads (≥Q20)
- Lower-quality but still polished consensus reads (<Q20)
- Unpolished consensus reads (`rq = -1`)
- 0- or 1-pass subreads unaltered (`rq = -1`)

If you want to only use HiFi reads, SMRT Link automatically generates additional
files for your convenience that only contain HiFi reads:

 - hifi_reads.**fastq**.gz
 - hifi_reads.**fasta**.gz
 - hifi_reads.**bam**

If you work with the `reads.bam` file directly, be aware that CCS reads of all
qualities are present. This file needs to be understood before piping
into your typical HiFi application.

## How does `--all` work?
With the special option `--all`, _ccs_ generates one representative
sequence per ZMW, irrespective of quality and passes.
For this, `--min-passes 0 --min-rq 0 --max-length 0` are set implicitly and
can't be changed; the maximum draft length filter is deactivated by this.
Filtering has to be performed downstream.

The workflow of _ccs_ with `--all` changes as follows.

**Exception 1:**
There is special behavior for low-pass ZMWs. If a ZMW has fewer than 2 full-length
subreads, use the subread of median length as representative consensus,
optionally with its kinetic information as forward orientation using `--all-kinetics`,
and do not polish.

**Exception 2:**
Only polish ZMWs with at least two full-length subreads mapping back to the draft.
Otherwise, set predicted accuracy `rq` tag to `-1` to indicate that the predicted
accuracy was not calculated and populate per-base QVs with `+` (QV10) the
approximate raw accuracy. Kinetic information are not available for unpolished
drafts.

**Exception 3:**
Instead of using an unpolished draft without kinetic information as representative
consensus sequence, if `--subread-fallback` is used, fall back to a
representative subread with kinetic information.

## How is `--all` different from explicitly setting `--min-passes 0 --min-rq 0`?
Setting `--min-passes 0 --min-rq 0` is the brute force combination that will
polish every ZMW, even those that only have one partial subread, with
polishing makes no difference.
In contrast, `--all` is a bit smarter and will only polish ZMWs with at
least one full-length subread and one additional partial subread; please see
previous paragraph how skipped ZMWs are represented in the output.

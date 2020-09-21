---
layout: default
parent: FAQ
title: Mode --all
---

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

### How is `--all` different from explicitly setting `--min-passes 0 --min-rq 0`?
Setting `--min-passes 0 --min-rq 0` is the brute force combination that will
polish every ZMW, even those that only have one partial subread, with
polishing makes no difference.
In contrast, `--all` is a bit smarter and will only polish ZMWs with at
least one full-length subread and one additional partial subread; please see
previous paragraph how skipped ZMWs are represented in the output.

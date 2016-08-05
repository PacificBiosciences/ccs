
``ConsensusCore2``: design and implementation
===========================================

Authors: David Alexander, Michael Brown, Nigel Delaney, Lance Hepler, Armin Töpfer
Last modified: August 5, 2016


Motivation
----------

PacBio raw reads are highly error prone compared to short read
technologies.  Confident inference about the underlying DNA sequence
can only be made on the basis of multiple observations---either
multiple read passes of the insert from a single ZMW, or multiple
reads otherwise determined to correspond to a common genomic region
(for example, because they map to the same window of the reference).
The problem of identifying the most plausible underlying DNA sequence,
based on multiple observations, is termed *consensus*.

Consensus is trivial for short read platforms because errors are
predominantly *mismatch errors*; in such a world one can confidently
call consensus by perfoming a multiple alignment and voting in each
column.

The consensus problem is more challenging for PacBio reads because the
predominant error type is insertions and deletions, especially in
homopolymer sequence.  While heuristic approaches based on multiple
alignment are possible even in this context, we have found better
results by taking the approach of modeling consenus in a probabilistic
framework, asking the question "what is the probability that an
underlying *template* sequence T generated an observed read sequence
R?", and then using the model :math:`\Pr(R \mid T)` to answer the
multi-read consensus problem via maximum likelihood---the consensus
sequence then mathematically defined as :math:`\arg\max_T \Pr(T \mid
\mathbf{R})`, where :math:`\mathbf{R}` represents the vector of
multiple reads.

Such a likelihood model is implemented using standard techniques from
the hidden Markov model (HMM) literature.  The likelihood model is
made tractable by approximating full dynamic programming using
banding, and, critically, the maximum likelihood search is made
tractable using a greedy search and a core routine using a
forward-backward identity that enables fast (:math:`O(1)`)
recalculation of the likelihood when the template is mutated pointwise
from :math:`T` to :math:`T+\mu`.

The consensus problem has multiple guises, as we have mentioned.  In
the *circular consensus sequence* application (*CCS*), multiple
successive reads from of a single circularized DNA insert template (a
*SMRTbell*) are gathered and used as evidence for a more confident
consensus sequence for the molecule.  In *genomic consensus*, reads
from many molecules are mapped to a reference genome or a "draft"
assembly, then a consensus sequence is generated for the genome by
aggregating together the reads mapping to each genomic region.  In the
context of a reference genome, this application is referred to as
*resequencing* and while a consensus sequence is produced, the goal is
typically to identify *variants*---positions evidencing variation from
the reference genome.  In the context of a "draft" assembly, the
application is referred to as "polishing" and the goal is simply to
produce the most accurate sequence of the genome as possible.

[add picture here]

Additionally, the mathematical model of :math:`Pr(T \mid \mathbf{R})`
has immediate applications even beyond the single-template consensus
problem.  It presents the natural framework for comparing hypotheses
about single or even mutliple underlying template sequences, giving
rise to elegant approaches for variant calling and confidence
estimation, as well as haplotype phasing.  The latter problem is


History
-------

An initial consensus model (names `Quiver` in 2012) was originated in
the original CCS implementation (ca. 2010) and proved capable of
generating consensus sequences with accuracy near Q25 provided enough
evidence---enough "passes" of the insert DNA.  The initial codebase
was in C# with some routines in C/C++.  In 2012-2013, the core of the
Quiver consensus algorithm was exported to a C++ library,
`ConsensusCore`, which provided SWIG bindings to higher-level host
languages (Python and C#, in chief), enabling use from other
applications .  The C# algorithm implementation was replaced by calls
to the `ConsensusCore` library; simultaneously, a Python application
(``GenomicConsenus``) was developed for applications in genome
assembly polishing, yielding very successful results: over Q50 (and in
some cases, over Q60) accuracy achieved on bacterial genome
assemblies.

The arrival of the ``GenomicConsensus" application coincided roughly
with the initial development of HGAP (the Hierarchical Genome Assembly
Program), thus offering fast end-to-end genome assembly with
exceptionally high quality results, establishing PacBio as a player in
the fields of "sequence finishing" and microbiology in 2013.
Throughput improvements would

Analyses of CCS results using the Quiver model showed evidence that
CCS was far from a solved problem.  First of all, the accuracy of CCS
results were found to fall short of the accuracy from multi-molecule
consensus sequences with comparable "coverage" in multiple controlled
experiments.  Perhaps more troubling, it was observed that consensus
accuracy "saturated", failing to increase beyond a certain point when


The Arrow model
---------------

Differences from Quiver...



Implementation
--------------

Draft consensus by partial-order alignment
``````````````````````````````````````````

The essence: calculating, and re-calculating, likelihood
````````````````````````````````````````````````````````

Numerical aspects
`````````````````

Scaling vs log-domain math.

The counterweights.


Setup for calling consensus: the `Integrator` classes
`````````````````````````````````````````````````````


Model lookup
````````````


Model parameter lookup
``````````````````````

Training
````````

Identifying (and removing) abberant reads: the "ZScore" concept
```````````````````````````````````````````````````````````````

Basic explanation goes here.

For the mathematical details of how the variance and expectation are
calculated, see :ref:`zscore-math`.



Appendices: the gory details
----------------------------

.. toctree::
   :maxdepth: 1

   ZScoreMath

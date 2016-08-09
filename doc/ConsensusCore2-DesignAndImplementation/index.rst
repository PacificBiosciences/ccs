
``ConsensusCore2``: design and implementation
===========================================


| **Authors:** David Alexander, Michael Brown, Nigel Delaney, Lance Hepler, Armin Töpfer
| **Last modified:** August 10, 2016


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
sequence then mathematically defined as :math:`\arg\max_T
\Pr(\mathbf{R} \mid T)`, where :math:`\mathbf{R}` represents the
vector of multiple reads.

Such a likelihood model is implemented using standard techniques from
the hidden Markov model (HMM) literature.  The likelihood model is
made tractable by approximating full dynamic programming using
banding, and, critically, the maximum likelihood search is made
tractable using a greedy search and a core routine using a
forward-backward identity that enables fast (:math:`O(1)`)
recalculation of the likelihood when the template is mutated pointwise,
:math:`T \to T+\mu`.

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

.. figure:: img/consensus-modalities.*
   :width: 100%

   *Consensus can applied in a single-molecule or multi-molecule
   context, for different applications.*

Additionally, the mathematical model of :math:`Pr(T \mid \mathbf{R})`
has immediate applications even beyond the single-template consensus
problem.  It presents the natural framework for comparing hypotheses
about single or even mutliple underlying template sequences, giving
rise to elegant approaches for variant calling and confidence
estimation, as well as haplotype phasing.  The latter problem is


History
-------

An initial consensus model (named `Quiver` in 2012) was originated in
the original CCS implementation (ca. 2010) and proved capable of
generating consensus sequences with accuracy near Q25 provided enough
evidence---enough "passes" of the insert DNA.  The initial codebase
was in C# with some routines in C/C++.  In 2012-2013, the core of the
Quiver consensus algorithm was exported to a C++ library,
``ConsensusCore``, which provided SWIG bindings to higher-level host
languages (Python and C#, in chief), enabling use from other
applications .  The C# algorithm implementation was replaced by calls
to the ``ConsensusCore`` library; simultaneously, a Python application
(``GenomicConsenus``) was developed for applications in genome
assembly polishing, yielding very successful results: over Q50 (and in
some cases, over Q60) accuracy achieved on bacterial genome
assemblies.

The arrival of the ``GenomicConsensus`` application coincided roughly
with the initial development of HGAP (the Hierarchical Genome Assembly
Program), thus offering fast end-to-end genome assembly with
exceptionally high quality results, establishing PacBio as a player in
the fields of "sequence finishing" and microbiology in 2013.
Throughput improvements over time then made PacBio assembly appealing
for larger genomes, including fungi, plants, and animals.

Analyses of CCS results using the Quiver model showed evidence that
CCS was far from a solved problem.  First of all, the accuracy of CCS
results were found to fall short of the accuracy from multi-molecule
consensus sequences with comparable "coverage" in multiple controlled
experiments.  Perhaps more troubling, it was observed that consensus
accuracy "saturated", failing to increase beyond a certain number of
passes.

Experiments suggested that part of this problem was due to a failure
to account for ZMW-specific variables in the Quiver model; "genomic"
consensus was not subject to this problem because multiple molecules
would counterbalance each other.  Another defect of the Quiver model
was that it was trained in a "discriminative" fashion, which we
believe biased the model in a manner that prevented convergence to
perfect accuracy as the number of CCS passes increased.

In addition to the deficiencies in accuracy, the Quiver model also was
burdensome from the a software engineering perspective.  It was
unusable as an "inference" tool---there was no way to use it to
estimate underlying physical sequencing HMM parameters (merge rates,
for example); a completely separate codebase and tool (EDNA) was used
for this purpose---imposing a maintenance burden.  Furthermore, the
nonstandard derivative-free optimization method that was used to train
Quiver was slow and unreliable, making it unsuitable for use in an
automated training pipeline.

These observations motivated the development of a new model, based on
standard likelihood theory and standard training techniques, and
suitable as an inference tool to replace EDNA.  This was the "UniteEM"
effort led by Nigel Delaney, resulting in the "Arrow" model.



The Arrow model
---------------

The Quiver model was effectively a "conditional random field" over
multivariate observations that encompassed the base calls and multiple
additional "QV" tracks emitted by the basecaller.  However it was not
trained in a standard manner; rather it was trained using an
unreliable and slow derivative-free optimizaton and with an objective
function that reflected a non-standard likelihood function more akin
to a classification accuracy.  Furthermore, its ability to adapt to
the

The Arrow model is a complete reboot.  It is a left-right HMM model,
very similar to the standard textbook left-right sequence "profile"
HMM.  Arrow differs from those models in a few key ways.  First, while
the standard profile HMM just models sequence alignment moves (Match,
Insert, Delete), Arrow models the enzymatic and photophysical events
(Sticks, Branches, (Mis)incorporations, Merges, and "dark" pulses)
underlying such moves in SMRT sequencing.  Secondly, emission and
transition parameters are not estimated independently for every state;
rather the states are "tied" by the dinucleotide template context.
Third, the transition parameters are "tied", depending only on the
template position and not on the state the model is in; this
simplification enables computing the model using a single matrix
matrix instead of one per state type.

.. figure:: img/hmms.*
   :width: 100%

   *The standard profile HMM (A) uses states to model alignment moves
   (Match, Insert, Delete), while the Arrow HMM (B) models the
   underlying physical events (which themselves give rise to
   "alignment" moves).  Otherwise, the models are quite similar.*

While the Arrow model can be used to estimate transition and emission
parameters freely for each read (the "EDNA" use case), in the more
common case emission parameters are only estimated per dinucleotide
context, combining all reads, while the transition probability
parameters are expected to follow a regression model on a few scalar
read covariates (at present, the only covariate used in this manner is
the read's SNR).

The Arrow model is trained using the standard Baum-Welch EM algorithm
procedure; the only extension here is that when the SNR is used as a
covariate, the M step requires maximum likelihood estimation of
a multinomial logistic model.

.. todo:: probably wise to include some actual math here!





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

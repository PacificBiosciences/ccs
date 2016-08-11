
.. _zscore-math:

"Z-Score": when is a read completely garbage?
------------------------

Garbage-In-Garbage-Out: It is not always useful to estimate fine-scale
distinctions when looking at examples that have excessively high
noise. Estimation is likely to be improved by simplying filtering away
all examples that are irrecoverably "broken" as long as you still have
sufficient sample size on which to make estimates.

With "single-molecule-weirdness" such as long bursty inserts,
individual reads (or subsections of them) might contain very high
levels of noise. The current record for bursty insert is about 5.7kb
of inserted sequence that has nothing to do with the reference
template being sequenced. This can throw off alignments and consensus
models (CCS2) that do not explicitly model these bursty insert
behaviours.

While it would be best to solve these behaviors upstream at the
chemical/polymerase level, we must have defenses in place to at least
identify when these undesirable behaviors present themselves so we
might at least try to filter them away before they wreak havoc on
estimates.

Z-Score for Filtering
---------------------

One attempt at identifying garbage reads is the "Z-Score" described
here [PBEP_4.pdf]

The idea is to compute the expected log probability of sequences
output by an HMM. Thus when presented with a sequence that has a
certain log probability, we can reason how far removed it is to
"normal" such that outliers can be filtered away.

Because an HMM is a Markov-independent series of transitions and
emissions each with corresponding probabilities, it is reasonable that
sums of log probabilities might be close to normal under certain
conditions.

Here is the basic HMM structure with match, branch (same match base
insert), stick (different match base insert), delete:

(pacbioHMM.svg.png)

Here is the result verbatim from PBEF_4.pdf:

| term       | meaning |
|------------| ------- |
| $p_m$      | match transistion probability |
| $p_d$      | delettion transistion probability |

# TODO: this doesn't seem to render correctly on github

LL = log(prob)

expected LL from context = expected LL from match/delete + expected LL from insertions

# break into simple match/delete and insert chains

expected LL from match/delete = $ (l_m + E[M]) \frac{p_m}{p_m+p_d} + l_d \frac{p_d}{p_m+p_d} $

# transistion weighted (LL from match trans+emis) and (LL from delete trans)

expected LL from insertions = $E[I] \frac{p_s+p_b}{p_m+p_d}$

# expected insertion LL weighted by expected length of insertion

$E[I] = (l_b+E[B]) \frac{p_b}{p_b+p_s} + (l_s+E[S]) \frac{p_s}{p_b+p_s}

# transition weighted LL branch transition/emission and LL stick transition/emission within insertion

As a sanity check, using a simple HMM, we generated random deviates,
computed means, and compared to the computed expected values.

The means agree nicely with a computed estimate of -101.7703 and an
empirical mean of -101.80 on 1000 random deviates for an HMM 120 bases
long. (The standard deviation is different 12.7 computed vs 16.01
empirial???)

Real-world performance on RSII data shows that the Z-score does have
good performance in filtering garbage reads.

Z-Score Shortcomings
--------------------

The bursty errors occur in localized regions. For a long read, these
localized bursts might not be detected by the Z-score metric.  Overall
the number of errors, if they were randomly distributed across the
read, might be within what might be expected normally. The fact that
they are all localized is what makes it abnormal.

An HMM can identify these localized bursts. The Viterbi path assigns
each match/delete state to a position in the read: (ref_i->read_j
prob_i). Because the HMM is a regular language, we known if ref_i
derives the string with prob_i and ref_{i+1} derives with prob_{i+1}
then ref_i derives it's portion with probability (prob_i / prob_{i+1}
or differences in log probability). This is the part of the HMM that
accounts for a single reference base. We can use the same Z-Score
ideas to determine outliers. If the subHMM derives 4 or less bases
99.999% of the time, then if in the Viterbi path, a derivation of 200
bases is observed, then we can conclude this is an outlier bursty
insert between this and then next reference base. Similar ideas exist
for forward / backward / posterior.

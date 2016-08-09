
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

Here is the result from PBEF_4.pdf:

| term       | meaning |
|------------| ------- |
| $p_m$      | match transistion probability |
| $p_d$      | delettion transistion probability |


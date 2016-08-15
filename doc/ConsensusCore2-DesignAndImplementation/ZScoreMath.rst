
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

The idea is to compute the expected mean and variance of the log
likelihood (LL) of sequences output by a given HMM. Thus when
presented with a sequence that has a certain log probability, we can
reason how far removed it is to "normal" such that outliers can be
filtered away.

An HMM is a series of Markov steps: transition and emisson. Each step
will additively change the LL. In a simple left-to-right profile HMM
with length, $T$, each match/delete will be visited once.  It is
reasonable that sums of log probabilities might be close to normal
under certain conditions.

Here is the basic HMM structure with match, branch (same match base
insert), stick (different match base insert), delete:

(pacbioHMM.svg.png)

For each "substate" (match, delete, branchIns, stickIns) we might
compute that substate's average contribtuion to the the LL over the
distribution of outputs induced by the HMM: E[LL] and Var[LL]. Then
$E[LL] = E[ \sum_{substates} LL] = \sum_{substates} E[LL]$ by
linearity of expectation. For IID substates this is $T*E[LL]$. Var[LL]
= $Var[ \sum_{substates} LL] = \sum_{substates} Var[LL]$ because each
substate is independent so all covariances are 0. For IID substates
this is $T*Var[LL]$.

For a motivational example, consider an "HMM" with a single
multinomial $p=(p_1,p_2,p_3,p_4)$. Here $E[LL] = p \dot
log(p). E[LL^2] = p \dot log^2(p). Var[LL]= E[LL^2]-E^2[LL]= p \dot
log^2(p) - (p \dot log(p))^2$. These are expected mean and variance of
the LL by this simple multinomial. By stringing together multinomial
emissions and transitions, we can estimate these quantities for an
HMM. Note the only non-trivial part is loops, but because we only
consider simple left-to-right profile HMMs, the only loops we
encounter are in the insert states.

Here is the result verbatim from PBEF_4.pdf:

| term       | meaning |
|------------| ------- |
| $p_m$      | match transistion probability |
| $p_d$      | deletion transistion probability |
| $p_s$      | stick transition probability |
| $p_b$      | branch transition probability |
| $l_m,l_d,l_s,l_b$ | log of transition probabilities |
| $E[\rightarrow]$ | expected LL from match or deletion |
| $E[\downarrow]$ | expected LL from all insertions | 
| $E[NN]$ | mean LL from context |
| $E[M]$ | mean logged match emission probability |
| $E[B]$ | mean logged branch emission probability |
| $E[S]$ | mean logged stick emission probability |
| $E[I]$ | mean logged transition,emission probability for insertion|

# TODO: this doesn't seem to render correctly on github

- $E[NN] = E[\downarrow] + E[\rightarrow]$

Break into simple match/delete and insert chains. Expected LL from
context = expected LL from match/delete + expected LL from insertions

- $E[\rightarrow] = $ (l_m + E[M]) \frac{p_m}{p_m+p_d} + l_d \frac{p_d}{p_m+p_d} $

Transistion weighted (LL from match trans+emis) and (LL from delete
trans)

- $E[\downarrow]= $E[I] \frac{p_s+p_b}{p_m+p_d}$

Expected insertion LL weighted by expected length of insertion where
$(p_s+p_b)$ is the probablity of looping in the insertion and
$(p_m+p_d)$ is the probability of looping out.

- $E[I] = (l_b+E[B]) \frac{p_b}{p_b+p_s} + (l_s+E[S]) \frac{p_s}{p_b+p_s}

Transition weighted LL branch transition/emission and LL stick
transition/emission within insertion

For the second moment, we can simply replace $LL^2$ for $LL$ in the
above equations everywhere except $E[\downarrow]$. For the
$E[\downarrow]$ insertions, the infinite sum of the geometric series
changes from the expectation.

Classic Identities:
| $E[X+Y]=E[X]+E[Y]$ |
| $Var(X) = E[X^2]-E^2[X]$ |
| $Var(X+Y) = Var(X) + Var(Y) + 2Cov(X,Y)$ |
| $E[a*X] = aE[X]$ |
| $Var(a*X) = a^2Var(X)$ |
| $E[XY]=E[X]E[Y]$ if independent |
| $\sum_{k=0}^\infty (1-p)^k*p*k*ll = ll*\frac{1-p}{p}$ |
| $\sum_{k=0}^\infty (1-p)^k*p*(k*ll)^2 = ll^2 \frac{(p-2)(p-1)}{p^2}$ |

Z-Score Sanity Check
--------------------

As a sanity check we generated random deviates using a simple HMM with
varying number of substates, computed means and variances, and
compared to the computed expected values.

The means are quite close computed versus estimated! However the
variance appears to be underestimated by the computation but still has
a linear relationship.

|  size |      mean |      var |  compmean |  compvar |
|-------|-----------|----------|-----------|----------|
|    32 | -26.65075 | 59.80562 | -27.13874 | 53.17624 |
|    60 | -50.99983 |  121.535 | -50.88515 | 99.70545 |
|   120 | -101.8052 | 256.4235 | -101.7703 | 199.4109 |
|   240 | -202.6475 | 473.0255 | -203.5406 | 398.8218 |
|   480 | -407.2655 | 1073.169 | -407.0812 | 797.6436 |

Real-world performance on RSII data shows that the Z-score does have
good performance in filtering garbage reads. Because we are able to
adjust the Z-Score threshold, good performance is obtained.

TODO: get the variance computation correct

Z-Score Shortcomings
--------------------

The bursty errors occur in localized regions. For a long read, these
localized bursts might not be detected by the Z-score metric.  Overall
the number of errors, if they were randomly distributed across the
read, might be within what might be expected normally. The fact that
they are all localized is what makes it abnormal.

An HMM can identify these localized bursts. The Viterbi path assigns
each match/delete state to a position in the read: $(ref_i->read_j
prob_i)$. Because the HMM is a regular language, we know if $ref_i$
derives the string starting at $read_j$ with $prob_i$ and $ref_{i+1}$
derives with $prob_{i+1}$ then $ref_i$ derives it's portion with
probability ($prob_i / prob_{i+1}$) or differences in log probability
if using log probability. This is the part of the HMM that accounts
for a single reference base. We can use the same Z-Score ideas to
determine outliers. If the substate HMM derives 4 or less bases
99.999% of the time, then if in the Viterbi path a derivation of 200
bases is observed, then we can conclude this is an outlier bursty
insert between this and then next reference base. (Similar ideas exist
for forward / backward / posterior.)

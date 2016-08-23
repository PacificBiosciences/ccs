
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


Z-Score Motivation
------------------

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
HMM.

An HMM emits a string by entering a state, emitting an output
according to an emission probability, making a random choice about
what state to transition to next, transitioning to that new state, and
repeating.

Consider an HMM state that has emission probabilities $e_j$,
transistion probabilities $t_i$, and new state $R_i$ that outputs the
rest of the output. The probability this state has is $(e_j t_i
\prod_k R_{ik})$ where the product is the rest of the states adding
their contributions. The log likelihood is $LL=\log e_j + \log t_i +
\sum_k \log R_{ik}$. Let $\log e = \sum_j e_j \log(e_j)$ and $\log^2 e
= \sum_j e_j \log^2(e_j)$ be the expectations over the possible
emission possibitilties.

The expectation of log likelihood is $EXP_{p_i} [(\log e + \log t_i +
\sum_k \log R_{ik})]$. Note the expecation of the log likelikhood of
the rest of the derivation represented by the sum is $EXP[ LL(R_i) ]$
so we get $EXP_{p_i} [(\log e + \log t_i + EXP[ LL(R_i)] ]$

The expectation of log likelihood squared is $EXP_{p_i}[ (\log e +
\log t_i + \sum_k \log R_{ik})^2]$. Expand the square. The only
complicating term is $\sum_k (\log R_{ik})^2 be we realize this to be
$EXP[ LL^2(R_i) ]$

Given this, we can write a dynamic program that computes the expected
log-likelihood and expected log-likelihood squared of $thisstate$
derving a string of $length$ for every state and every length up to
some maximum.

The relvant python code that computes the expectation over the
different choices of next rest state (nexts) with transition
probabilities (nextp) and expected emissions (expEmLL, expEmLL2):

'''python
        # cycle through choices
        overallresult = [0.0, 0.0, 0.0]
        for ch in range(len(this.nexts)):

            transp = this.nextp[ch]
            rhs = self.derive( this.nexts[ch], newlength)

            A = self.mylog(transp) # log transition
            A2 = A*A               # log transistion squared
            B = expEmLL            # expected log emission likelihood
            B2= expEmLL2           # expected log squared emission likelihood
            C = rhs[1]             # EXP[LL] of next
            C2= rhs[2]             # EXP[LL^2] of next
            thisll =  transp*(A+B+C)
            thisll2 = transp*( A2 + B2 + C2 + 2*A*B + 2*A*C +2*B*C)

            thisprob = transp*rhs[0]

            if thisprob>0.0:
                overallresult[0]=overallresult[0]+thisprob
                overallresult[1]=overallresult[1]+thisll
                overallresult[2]=overallresult[2]+thisll2

        self.memo[self.key(thisstate,length)]=overallresult
        return(overallresult)
'''

Z-Score for Filtering
---------------------

One method for identifying garbage reads is the "Z-Score" described
here [PBEP_4.pdf]. This result sums out the infinite sums caused by
the looping insert states analytically to get an analytic result
rather than a dynamic programming one.

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

- $E[I] = (l_b+E[B]+E[\rightarrow]) \frac{p_b}{p_b+p_s} + (l_s+E[S]+E[\rightarrow]) \frac{p_s}{p_b+p_s}

Transition weighted LL branch transition/emission and LL stick
transition/emission within insertion. (Note this updates to three
terms versus to two in the PBEP)

For the second moment, we can simply replace $LL^2$ for $LL$ in the
above equations.

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

The means and variances are close computed versus estimated.

|  size |     mean |     var | compmean | compvar|
|-------|----------|---------|----------|--------|
|    32 |-26.65075 |59.80562 |-27.13874 |68.15969|
|    60 |-50.99983 | 121.535 |-50.88515 |127.7994|
|   120 |-101.8052 |256.4235 |-101.7703 |255.5988|
|   240 |-202.6475 |473.0255 |-203.5406 |511.1977|
|   480 |-407.2655 |1073.169 |-407.0812 |1022.395|


Real-world performance on RSII data shows that the Z-score does have
good performance in filtering garbage reads. Because we are able to
adjust the Z-Score threshold, good performance is obtained.

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

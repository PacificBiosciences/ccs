<p align="center">
  <img src="img/ccs.png" alt="CCS logo" width="250px"/>
</p>
<h1 align="center">CCS</h1>
<p align="center">Generate Highly Accurate Single-Molecule Consensus Reads (HiFi Reads)</p>

***

_ccs_ takes multiple subreads of the same SMRTbell molecule and combines
them using a statistical model to produce one highly accurate consensus sequence,
also called HiFi read, with base quality values.
This tool powers the _Circular Consensus Sequencing_ workflow in SMRT Link.

## Availability
Latest `ccs` can be installed via bioconda package `pbccs`.

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
for information on Installation, Support, License, Copyright, and Disclaimer.

## Latest Version
Version **4.2.0**: [Full changelog here](#changelog)

## Schematic Workflow
<p align="center"><img width="600px" src="img/ccs-workflow.png"/></p>

## Execution
**Input**: Subreads from a single movie in PacBio BAM format (`.subreads.bam`).

**Output**: Consensus reads in a format inferred from the file extension:
unaligned BAM (`.bam`); bgzipped FASTQ (`.fastq.gz`);
or SMRT Link XML (`.consensusreadset.xml`) which also generates a corresponding
BAM file.

Run on a full movie:

    ccs movie.subreads.bam movie.ccs.bam

Parallelize by using `--chunk`.
See [how-to chunk](#how-can-I-parallelize-on-multiple-servers).

## FAQ

### What impacts the number and quality of HiFi reads that are generated?
The longer the polymerase read gets, more passes of the SMRTbell
are produced and consequently more evidence is accumulated per molecule.
This increase in evidence translates into higher consensus accuracy, as
depicted in the following plot:

<p align="center"><img width="600px" src="img/ccs-acc.png"/></p>

### How is number of passes computed?
Each read is annotated with a `np` tag that contains the number of
full-length subreads used for polishing. Full-length subreads are flanked by
adapters and thus cover the full insert.
Since the first version of _ccs_, number of passes has only accounted for
full-length subreads. In version v3.3.0 windowing has been added, which
takes the minimum number of full-length subreads across all windows.
Starting with version v4.0.0, minimum has been replaced with mode to get a
better representation across all windows. Only subreads that pass the subread
length filter (please see next FAQ about filters) and were not dropped during
polishing are accounted for.

Similarly, tag `ec` reports effective coverage, the average subread coverage
across all windows. This metric includes all subreads, independent of being
full- or partial-length subreads, that pass length filters and did not fail
during polishing. In most cases `ec` will be roughly `np + 1`.

### Which and in what order are filters applied?
_ccs_ exposes the following filters on input subreads, draft consensus,
and final output consensus:

    Input Filter Options:
      --min-passes   INT    Minimum number of full-length subreads required to generate CCS for a ZMW. [3]
      --min-snr      FLOAT  Minimum SNR of subreads to use for generating CCS [2.5]

    Draft Filter Options:
      --min-length   INT    Minimum draft length before polishing. [10]
      --max-length   INT    Maximum draft length before polishing. [50000]

    Output Filter Options:
      --min-rq       FLOAT  Minimum predicted accuracy in [0, 1]. [0.99]

Data flow how each ZMW gets processed and filtered:
1. Remove subreads with lengths <50% or >200% of the median subread length.
2. Remove subreads with SNR below `--min-snr`.
3. Stop if number of full-length subreads is fewer than `--min-passes`.
4. Generate draft sequence and stop if draft length does not pass `--min-length` and `--max-length`.
5. Polish consensus sequence and only emit read if predicted accuracy is at least `--min-rq`.

**Regarding the subread length filter, why did we pick 50% and 200%?**\
If you sequence a palindromic sequence and one of your subreads is smaller than
50%, then it becomes a challenge to determine where exactly the subread fits
without alignment hacks; similarly a short subread might map wrongly to a highly
repetitive sequence.\
Subreads that are longer than 200% the median length are likely due to missed
adapter calls, where the initial polymerase read was wrongly split into subreads.

### How do I read the ccs_report.txt file?
The `ccs_report.txt` file summarizes (B) how many ZMWs generated HiFi reads (CCS) and
(C) how many failed CCS generation because of the listed causes. For (C), each ZMW
contributes to exactly one reason of failure; percentages are with respect to (C).

The following comments refer to the filters that are explained in the FAQ above.

    ZMWs input          (A)  : 44327
    ZMWs generating CCS (B)  : 16304 (36.78%)
    ZMWs filtered       (C)  : 28023 (63.22%)

    Exclusive ZMW counts for (C):
    Median length filter     : 0 (0.00%)      <- All subreads were filtered in (1)
    Below SNR threshold      : 574 (2.05%)    <- All subreads were filtered in (2)
    Lacking full passes      : 24706 (88.16%) <- Fewer than --min-passes full-length (FL) reads (3)
    Heteroduplex insertions  : 422 (1.51%)    <- Single-strand artifacts
    Coverage drops           : 29 (0.10%)     <- Coverage drops would lead to unreliable polishing results
    Insufficient draft cov   : 282 (1.01%)    <- Not enough subreads aligned to draft end-to-end
    Draft too different      : 0 (0.00%)      <- Fewer than --min-passes FL reads aligned to draft
    Draft generation error   : 94 (0.34%)     <- Subreads don't agree to generate a draft sequence
    Draft above --max-length : 21 (0.07%)     <- Draft sequence is longer than --min-length (4)
    Draft below --min-length : 0 (0.00%)      <- Draft sequence is shorter than --min-length (4)
    Reads failed polishing   : 0 (0.00%)      <- Too many subreads were dropped while polishing
    Empty coverage windows   : 0 (0.00%)      <- At least one window has no coverage
    CCS did not converge     : 0 (0.00%)      <- Draft has too many errors that can't be polished in time
    CCS below minimum RQ     : 1916 (6.84%)   <- Predicted accuracy is below --min-rq (5)
    Unknown error            : 0 (0.00%)      <- Rare implementation errors

If run in `--by-strand` mode, rows may contain half ZMWs, as we account
each strand as half a ZMW.

### What is the definition of a heteroduplex?
In general, whenever bases on one strand of the SMRTbell are not the
reverse complement of the other strand, as small as a single base `A` with a
matching `G`. _ccs_ would polish this to one of the bases and reflect the
ambiguity in the base QV. In our case, when one strand has more than `20`
additional bases that the other strand does not have, _ccs_ won't be able to
converge to a consensus sequence, and consequently will remove the ZMW and
increase the counter for heteroduplexes found in the `ccs_report.txt` file.

### How can I parallelize on multiple servers?
You can parallelize by chunking. Since _ccs_ v4.0.0, direct chunking via `--chunk`
is possible. For this, the `.subreads.bam` file must accompanied by a
`.pbi` file. To generate the index `subreads.bam.pbi`, use
`pbindex`, which can be installed with `conda install pbbam`.

    pbindex movie.subreads.bam

An example workflow, all ccs invocations can run simultaneously:

    ccs movie.subreads.bam movie.ccs.1.bam --chunk 1/10 -j <THREADS>
    ccs movie.subreads.bam movie.ccs.2.bam --chunk 2/10 -j <THREADS>
    ...
    ccs movie.subreads.bam movie.ccs.10.bam --chunk 10/10 -j <THREADS>

Merge chunks with `pbmerge` and index with `pbindex`

    pbmerge -o movie.ccs.bam movie.ccs.*.bam
    pbindex movie.ccs.bam

or use `samtools`

    samtools merge -@8 movie.ccs.bam movie.ccs.*.bam

### What happened to unanimity?
Unanimity lives on as a PacBio internal library to generate consensus sequences.
Documentation will be limited to _ccs_ distributed via bioconda.

### Where is the source code?
We have stopped mirroring code changes to GitHub in March 2018.
Instead, we provide binaries on bioconda to ease end-user experience.
If your project relies on outdated unanimity source code,
please use [this commit](https://github.com/PacificBiosciences/unanimity/tree/6f11a13e1472b8c00337ba8c5e94bf83bdab31d6).

### Help! I am getting "Unsupported ..."!
If you encounter the error `Unsupported chemistries found: (...)` or
`unsupported sequencing chemistry combination`, your _ccs_ binaries do not
support the used sequencing chemistry kit, from here on referred to as "chemistry".
This may be because we removed support of an older chemistry or your binary predates
release of the used chemistry.
This is unlikely to happen with _ccs_ from SMRT Link installations, as SMRT Link
is able to automatically update and install new chemistries.
Thus, the easiest solution is to always use _ccs_ from the SMRT Link version that
shipped with the release of the sequencing chemistry kit.

**Old chemistries:**
With _ccs_ 4.0.0, we have removed support for the last RSII chemistry `P6-C4`.
The only option is to downgrade _ccs_ with `conda install pbccs==3.4`.

**New chemistries:**
It might happen that your _ccs_ version predates the sequencing chemistry kit.
To fix this, install the latest version of _ccs_ with `conda update --all`.
If you are an early access user, follow the [monkey patch tutorial](#monkey-patch-ccs-to-support-additional-sequencing-chemistry-kits).

### Monkey patch _ccs_ to support additional sequencing chemistry kits
Please create a directory that is used to inject new chemistry information
into _ccs_:

```sh
mkdir -p /path/to/persistent/dir/
cd /path/to/persistent/dir/
export SMRT_CHEMISTRY_BUNDLE_DIR="${PWD}"
mkdir -p arrow
```

Execute the following step by step instructions to fix the error you are observing
and afterwards proceed using _ccs_ as you would normally do. Additional chemistry
information is automatically loaded from the `${SMRT_CHEMISTRY_BUNDLE_DIR}`
environmental variable.

#### Error: "unsupported sequencing chemistry combination"
Please download the latest out-of-band `chemistry.xml`:

```sh
wget https://raw.githubusercontent.com/PacificBiosciences/pbcore/develop/pbcore/chemistry/resources/mapping.xml -O "${SMRT_CHEMISTRY_BUNDLE_DIR}"/chemistry.xml
```

#### Error: "Unsupported chemistries found: (...)"
Please get the latest consensus model `.json` from PacBio and
copy it to:

```sh
cp /some/download/dir/model.json "${SMRT_CHEMISTRY_BUNDLE_DIR}"/arrow/
```

### How fast is CCS?
We tested CCS runtime using 500 ZMWs per length bin with exactly 7 passes.

<img width="1000px" src="img/runtime.png"/>


#### How does that translate into time to result per SMRT Cell?
We will measure time to result for Sequel I and II CCS sequencing collections
on a PacBio recommended HPC, according to the
[Sequel II System Compute Requirements](https://www.pacb.com/wp-content/uploads/SMRT_Link_Installation_v701.pdf)
with 192 physical or 384 hyper-threaded cores.

1) Sequel System: 15 kb insert size, 30-hours movie, 37 GB raw yield, 2.3 GB HiFi UMY
2) Sequel II System: 15 kb insert size, 30-hours movie, 340 GB raw yield, 24 GB HiFi UMY

CCS version | Sequel System | Sequel II System
:-: | :-: | :-:
≤3.0.0 | 1 day | >1 week
3.4.1 | 3 hours | >1 day
4.0.0 | 40 minutes | 6 hours
≥4.2.0 | **30 minutes** | **4 hours**

#### How is CCS speed affected by raw base yield?
Raw base yield is the sum of all polymerase read lengths.
A polymerase read consists of all subreads concatenated
with SMRTbell adapters in between.

Raw base yield can be increased with
1) higher percentage of single loaded ZMWs and
2) longer movie times that lead to longer polymerase read lengths.

Since the first version, _ccs_ scaled linear in (1) the number of single loaded
ZMWs per SMRT Cell.
Starting with version 3.3.0 _ccs_ scaled linear in (2) the polymerase read length
and with version 4.0.0 _ccs_ scales sublinear.

#### What did change in each version?
CCS version | O(insert size) | O(#passes)
:-: | :-: | :-:
≤3.0.0 | quadratic | linear
3.4.1 | **linear** | linear
≥4.0.0 | linear | **sublinear**

#### How can version 4.0.0 be sublinear in the number of passes?
With the introduction of new heuristics, individual draft bases can skip
polishing if they are of sufficient quality.
The more passes a ZMW has, the fewer bases need additional polishing.

### What heuristics are used?
Following heuristics are enabled
 - determine which bases need polishing,
 - remove ZMWs with single-strand artifacts such as heteroduplexes
 - remove large insertions that likely are due to sequencing errors,
 - on-the-fly model caching with less SNR resolution,
 - adaptive windowing strategy with a target window size of 22 bp with ±2 bp overlaps, avoiding breaks in simple repeats (homopolymers to 4mer repeats)

### Does speed impact quality and yield?
Yes it does. With ~35x speed improvements from version 3.1.0 to 4.0.0 and
consequently reducing CPU time from >60,000 to <2,000 core hours,
heuristics and changes in algorithms lead to slightly lower yield and
accuracy if run head to head on the same data set. Internal tests show
that _ccs_ 4.0.0 introduces no regressions in CCS-only Structural Variant
calling and has minimal impact on SNV and indel calling in DeepVariant.
In contrast, lower DNA quality has a bigger impact on quality and yield.

### Can I tune _ccs_ to get improved results?
No, we optimized _ccs_ such that there is a good balance between speed and
output quality.

### Can I produce one consensus sequence for each strand of a molecule?
Yes, please use `--by-strand`. Make sure that you have sufficient coverage,
as `--min-passes` are per strand in this case. For each strand, _ccs_
generates one consensus read that has to pass all filters.
Read name suffix indicates strand. Example:

    m64011_190714_120746/14/ccs/rev
    m64011_190714_120746/35/ccs/fwd

How does `--by-strand` work? For each ZMW:
 * Determine orientation of reads with respect to the one closest to the mean length
 * Sort reads into two buckets, forward and reverse strands
 * Treat each strand as an individual entity as we do with ZMWs
   * Apply all filters per strand individually
   * Create a draft for each strand
   * Polish each strand
   * Write out each polished strand consensus

### Is there a progress report?
Yes. With `--log-level INFO`, _ccs_ provides status to `stderr` every
`--refresh-rate seconds` (default 30):

    70/420/26.3 19/114/7.1 13s

The log output explains each field:

    Logging info: Z1/Z2/Z3 C1/C2/C3 ETA [X]
    Z1: #ZMWs processed since start
    Z2: #ZMWs processed in the last minute
    Z3: #ZMWs processed in the last minute per thread
    C1: #CCSs generated since start
    C2: #CCSs generated in the last minute
    C3: #CCSs generated in the last minute per thread
    ETA: Estimated remaining processing time
    [X]: Per minute metrics are extrapolated

If there is no `.pbi` file present, ETA will be omitted.

### The binary does not work on my linux system!
Contrary to official SMRT Link releases, binaries distributed via bioconda are
tuned for performance while sacrificing backward compatibility.
We are aware of following errors and limitations. If yours is not listed, please
file an issue on our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda).

**`Illegal instruction`** Your CPU is not supported.
A modern (post-2008) CPU with support for
[SSE4.1 instructions](https://en.wikipedia.org/wiki/SSE4#SSE4.1) is required.
SMRT Link also has this requirement.

**`FATAL: kernel too old`** Your OS or rather your kernel version is not supported.
Since CCS v4.2 we also ship a second binary via bioconda `ccs-alt`, which uses
`mimalloc` as a modern allocator instead of bundling a newer `glibc`. Please use
this alternative binary.

### Why do I get more yield if I increase `--min-passes`
For versions newer than 3.0.0 and older than 4.2.0, we required that after
draft generation, at least `--min-passes` subreads map back to the draft.
Imagine the following scenario, a ZMW with 10 subreads generates a draft to which
only a single subread aligns. This draft is of low quality and does not
represent the ZMW, yet if you ask for `--min-passes 1`, this low quality draft
is being used. Starting with version 4.2.0, we switch to an additional
percentage threshold of 51% aligning subreads to avoid this problem. This fixes
the majority of discrepencies for less than three passes.

Why do we have this problem at all, shouldn't the draft stage be robust enough?
Robustness comes with speed trade-offs. We have a cascade of different draft
generators, from very fast and unstabler to slow and robust. If a ZMW fails
to generate a draft for a fast generator, it falls back multiple times until it
reaches the slower and more robust generator. This approach is still much faster
than always relying on the robust generator.

## Licenses
PacBio® tool _ccs_, distributed via Bioconda, is licensed under
[BSD-3-Clause-Clear](https://spdx.org/licenses/BSD-3-Clause-Clear.html)
and statically links GNU C Library v2.29 licensed under [LGPL](https://spdx.org/licenses/LGPL-2.1-only.html).
Per LPGL 2.1 subsection 6c, you are entitled to request the complete
machine-readable work that uses glibc in object code.


## Changelog

 * 4.2.0:
   * SMRT Link v9.0 release
   * Speed improvements
   * Minor yield improvements, by requiring a percentage of subreads mapping back to draft instead of `--min-passes`
   * Add effective coverage `ec` tag
   * Lowering `--min-passes` does no longer reduce yield
   * Add `--batch-size` to better saturate machine with high core counts
   * Simplify log output
   * Fix bug in predicted accuracy calculation
   * Improved `ccs_report.txt` summary
 * 4.1.0:
   * Minor speed improvements
   * Fix `--by-strand` logic, see more [here](#can-i-produce-one-consensus-sequence-for-each-strand-of-a-molecule)
   * Allow vanilla `.xml` output without specifying dataset type
   * Compute wall start/end for each output read (future basecaller functionality)
 * 4.0.0:
   * SMRT Link v8.0 release
   * Speed improvements
   * Removed support for legacy python Genomic Consensus, please use [gcpp](https://github.com/PacificBiosciences/gcpp)
   * New command-line interface
   * New report file
 * 3.4.1
   * Released with SMRT Link v7.0
   * Log used chemistry model to INFO level
 * 3.4.0
   * Fixes to unpolished mode for IsoSeq
   * Improve runtime when `--minPredictedAccuracy` has been increased
 * 3.3.0
   * Add a windowing approach to reduce computational complexity from quadratic to linear
   * Improve multi-threading framework to increase throughput
   * Enhance XML output, propagate `CollectionMetadata`
   * Includes latest chemistry parameters
 * 3.1.0
   * Add `--maxPoaCoverage` to decrease runtime for unpolished output, special parameter for IsoSeq workflow
   * Chemistry parameters for SMRT Link v6.0

## DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

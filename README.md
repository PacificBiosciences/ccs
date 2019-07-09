<p align="center">
  <img src="img/ccs.png" alt="CCS logo" width="250px"/>
</p>
<h1 align="center">CCS</h1>
<p align="center">Generate <span style=font-weight:bold>Hi</span>gh-<span style=font-weight:bold>Fi</span>delity Single-Molecule Consensus Reads</p>

***

_ccs_ takes multiple (sub)reads of the same SMRTbell molecule and combines
them using a statistical model to produce one high-fidelity consensus sequence
with base quality values.
This tool powers the _Circular Consensus Sequences_ workflow in SMRT Link.

## Availability
Latest `ccs` can be installed via bioconda package `pbccs`.

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
for information on Installation, Support, License, Copyright, and Disclaimer.

## Latest Version
Version **3.4.1**: [Full changelog here](#changelog)

## Schematic Workflow
<p align="center"><img width="600px" src="img/ccs-workflow.png"/></p>

## Execution
**Input**: Subreads from a single movie in PacBio BAM format (`.subreads.bam`).

**Output**: Consensus reads in a format inferred from the file extension:
unaligned BAM (`.bam`); FASTQ (`.fastq`);
or SMRT Link XML (`.consensusreadset.xml`) which also generates a corresponding
BAM file.

Run on a full movie:

    ccs movie.subreads.bam movie.ccs.bam

## FAQ

### What impacts the number and quality of CCS reads that are generated?
The longer the polymerase read gets, more readouts (passes) of the SMRTbell
are produced and consequently more evidence is accumulated per molecule.
This increase in evidence translates into higher consensus accuracy, as
depicted in the following sketch:

<p align="center"><img width="600px" src="img/ccs-acc.png"/></p>

### What happened to unanimity?
Unanimity lives on as a PacBio internal library to generate consensus sequences.
Customer-facing documentation will be limited to _ccs_ distributed via bioconda.

### Where is the source code?
We have stopped mirroring code changes to GitHub in March 2018.
Instead, we provide binaries on bioconda to ease end-user experience.
If your project relies on outdated unanimity source code,
please use [this commit](https://github.com/PacificBiosciences/ccs/tree/6f11a13e1472b8c00337ba8c5e94bf83bdab31d6).

### Help! I am getting "Unsupported ..."!
If you encounter the error `Unsupported chemistries found: (...)` or
`unsupported sequencing chemistry combination`, your _ccs_ binaries predates
the used sequencing chemistry kit.
This is unlikely to happen with _ccs_ from SMRT Link installations, as SMRT Link
is able to automatically update and install new chemistries.

If you are not an early access user, stop here and solve it via updating your
_ccs_ with `conda update --all`.

If you are an early access user, please create a directory that is used to
inject new chemistry information into _ccs_:

```sh
mkdir -p /path/to/persistent/dir/
cd /path/to/persistent/dir/
export SMRT_CHEMISTRY_BUNDLE_DIR="${PWD}"
mkdir -p arrow
```

To fix `unsupported sequencing chemistry combination` please download the latest
out-of-band `chemistry.xml`:

```sh
wget https://raw.githubusercontent.com/PacificBiosciences/pbcore/develop/pbcore/chemistry/resources/mapping.xml -O "${SMRT_CHEMISTRY_BUNDLE_DIR}"/chemistry.xml
```

To fix `Unsupported chemistries found: (...)` please get the latest consensus
model `.json` from PacBio and copy it to:

```sh
cp /some/download/dir/model.json "${SMRT_CHEMISTRY_BUNDLE_DIR}"/arrow/
```

Afterwards, proceed using _ccs_ as you would normally do. New chemistry
information is automatically loaded from the `${SMRT_CHEMISTRY_BUNDLE_DIR}`
environmental variable.

## Changelog

 * **3.4.1**
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

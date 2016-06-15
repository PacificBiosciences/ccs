<h1 align="center">
    pbccs - Generate Accurate Consensus Sequences from a Single SMRTbell
    <a href="https://circleci.com/gh/PacificBiosciences/pbccs">
        <img src="https://circleci.com/gh/PacificBiosciences/pbccs.svg?style=svg"
             alt="Circle CI" />
    </a>
</h1>
<p align="center">
  <img src="http://www.evolvedmicrobe.com/CCS.png" alt="Image of SMRTbell"/>
</p>

The ccs program takes multiple reads of the same SMRTbell sequence and combines 
them, employing a statistical model, to produce one high quality consensus sequence.

## Overview

 - [Installation](BUILD.md)
 - [Usage](USAGE.md)
  - [Input](USAGE.md#input)
  - [Running CCS](USAGE.md#running-ccs)
  - [Output](USAGE.md#output)
  - [QUAL values](USAGE.md#interpretting-qual-values)
  - [Z-Scores](USAGE.md#what-are-z-scores)
  - [CCS Yield Report](USAGE.md#ccs-yield-report)

## Quick Start

Convert legacy RSII h5 data to BAM

    bax2bam -o movieName movieName.1.bax.h5 movieName.2.bax.h5 movieName.3.bax.h5

Compute CCS reads from Sequel or RSII BAM file

    ccs myresults.bam movieName.subreads.bam

## Help

Futher help can be requested with -h. Issues? Bugs? Please create a github issue.
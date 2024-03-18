---
layout: default
title: Changelog
nav_order: 99
---

# Version changelog

**8.0.1**
   * SMRT Link v13.1 release
   * No customer-facing changes

8.0.0
   * SMRT Link v13.0 release
   * Always encode HiFi kinetics with v1 codec
   * [Add fail flag `ff`](/faq/fail-reads)

7.0.0
   * SMRT Link v12.0 release
   * Includes DeepConsensus polishing on the Revio platform
   * Arrow polishing on GPU on the Revio platform
   * New output files layout
   * [Output BAM file compression](/faq/qv-binning)

6.4.0
   * [Single-strand HiFi kinetics](/faq/kinetics)
   * Faster draft generation
   * CLR subsampling stores XML and PBI

6.3.0
   * SMRT Link v11.0 release
   * Heteroduplex finder
      * Change file output, single file with multiple read groups
      * Add HD detection for substitutions differences between strands
   * Add missing adapter tags `ma` and `ac`
   * Adhere to latest [PacBio BAM spec](https://pacbiofileformats.readthedocs.io/)
   * Add `--subsample-clr-perc` and `--subsample-clr-file` to store a percentage of productive ZMWs as subreads
   * Add `--fastq` as an additional output file

6.2.0
   * SMRT Link v10.2 release
   * Improved low-complexity handling, better runtime and lower memory usage
   * Improved BAM merge step, up to 5x faster
   * Improved compute run time
   * `INFO` logging if chemistry bundle is injected
   * Enable [strand splitting](/faq/mode-heteroduplex-filtering) of ZMWs that contain large insertion heteroduplexes
   * Use `TMPDIR` environment variable for storing temporary files
   * New `INFO` log summary and `ccs_reports.txt`

6.0.0
   * SMRT Link v10.1 release
   * Increase number of HiFi reads
   * Increase percentage of barcode yield
   * Run time, CPU time, and peak RSS improvements
   * Change main draft algorithm from pbdagcon to sparc
   * Replace minimap2 with pancake and edlib/KSW2

5.0.0
   * SMRT Link v10.0 release
   * Add `--hifi-kinetics` to average kinetic information for polished reads
   * Add `--all-kinetics` to add kinetic information for all ZMWs, except for unpolished draft consensus
   * Add `--subread-fallback`, combined with `--all`, use a subread instead of a draft as representative consensus
   * Use sDUST to identify tandem repeats
   * Output HiFi yield (>= Q20) and Unique Molecular Yield as INFO log
   * Set `--top-passes 60` default
   * Abort if chemistry information is missing in BAM header
   * Add non-blocking temporary file writing
   * Add `--input-buffer` to smooth IO fluctations
   * Add `--all` to generate one representative read per ZMW
   * Reuse prefix of output file for report files to avoid unintentional clobbering
   * Add `zmw_metrics.json`, metrics about each ZMW; file name can be set with `--metrics-json`
   * Add JSON output of ccs_reports via `--report-json`
   * Add `--suppress-reports` to suppress generating default report and metric files

4.2.0
   * SMRT Link v9.0 release
   * Speed improvements
   * Minor yield improvements, by requiring a percentage of subreads mapping back to draft instead of `--min-passes`
   * Add effective coverage `ec` tag
   * Lowering `--min-passes` does no longer reduce yield
   * Add `--batch-size` to better saturate machine with high core counts
   * Simplify log output
   * Fix bug in predicted accuracy calculation
   * Improved `ccs_report.txt` summary

4.1.0
   * Minor speed improvements
   * Fix `--by-strand` logic, see more [here](https://ccs.how/faq/mode-by-strand)
   * Allow vanilla `.xml` output without specifying dataset type
   * Compute wall start/end for each output read (future basecaller functionality)

4.0.0
   * SMRT Link v8.0 release
   * Speed improvements
   * Removed support for legacy python Genomic Consensus, please use [gcpp](https://github.com/PacificBiosciences/gcpp)
   * New command-line interface
   * New report file

3.4.1
   * SMRT Link v7.0 release
   * Log used chemistry model to INFO level

3.4.0
   * Fixes to unpolished mode for IsoSeq
   * Improve runtime when `--minPredictedAccuracy` has been increased

3.3.0
   * Add a windowing approach to reduce computational complexity from quadratic to linear
   * Improve multi-threading framework to increase throughput
   * Enhance XML output, propagate `CollectionMetadata`
   * Includes latest chemistry parameters

3.1.0
   * Add `--maxPoaCoverage` to decrease runtime for unpolished output, special parameter for IsoSeq workflow
   * Chemistry parameters for SMRT Link v6.0

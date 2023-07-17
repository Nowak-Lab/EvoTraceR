EvoTraceR
================

[![Actions Status](https://github.com/Nowak-Lab/EvoTraceR/workflows/check-master/badge.svg)](https://github.com/Nowak-Lab/EvoTraceR/actions?query=workflow%3Acheck-master)
[![Actions Status](https://github.com/Nowak-Lab/EvoTraceR/workflows/check-development/badge.svg)](https://github.com/Nowak-Lab/EvoTraceR/actions?query=workflow%3Acheck-development)

*EvoTraceR* is an R package to analyse sequencing amplicon data from CRISPR-Cas9 recorder lineage tracing experiments. The package takes in paired-end FASTQ files from one to many tissues. The sequenced amplicon can contain one to many Cas9 cut sites. *EvoTraceR* trims and merges reads, collapses duplicates, calls mutations and infers a tree using [Cassiopeia](https://github.com/YosefLab/Cassiopeia). The package outputs the inferred tree of relationships between Amplicon Sequence Variants (ASVs) as well as summary plots and tables of mutations.

Please feel free to contact us with feedback: Dawid Nowak, dgn2001 at med.cornell.edu or Armin Scheben, ascheben at cshl.edu.

Installation
--------------

```
library(devtools)
install_github("Nowak-Lab/EvoTraceR")
```




# AACR Project GENIE

## Introduction

This repository documents code used to gather, QC, standardize, and analyze data uploaded by institutes participating in AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange). Follow instructions below to download and populate the folder `data` with the core Genie dataset. And then review the contents of the `analyses` folder, which contains a sub-folder each, for various analyses that can be performed on those datasets.

## Dependencies

These are tools or packages you will need, to be able to reproduce these results:
- A Linux machine or a cluster of compute nodes with a job scheduler like LSF (`bsub`) or SGE (`qsub`)
- Python 2.7.10 or higher
- Sage Synapse's [command-line client](http://python-docs.synapse.org/CommandLineClient.html) (`pip install synapseclient`)

## Fetch Data

Download the core dataset into `data/cbioportal`. This data feeds the web-based UI at [cbioportal.org/genie](http://www.cbioportal.org/genie/):
```bash
synapse get -r --downloadLocation data/cbioportal syn5521835
```

## Reporting bugs

Please double check and triple check your bug, and if you're super certain that we're not infallible, then [click here](https://github.com/Sage-Bionetworks/Genie/issues) to report the bug, along with a clear explanation of how we can find or reproduce it.

# AACR Project GENIE

## Introduction

This repo documents code used to gather, QC, and standardize data uploaded by participating sites for AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange). Follow instructions below to create and populate a folder for `data`. And then review the contents of the `analyses` folder, which contains a subfolder each, for various analyses that could be performed on those datasets.

## Dependencies

These are tools or packages you will need, to be able to reproduce these results:
- A Linux machine or a cluster of compute nodes with a job scheduler like `LSF` or `SGE`
- Python 2.7.10 or higher
- Sage Synapse's [command-line client](http://python-docs.synapse.org/CommandLineClient.html) (`pip install synapseclient`)

## Fetch Data

Download cBioPortal data into a subfolder named `data/cbioportal`:
```bash
synapse -u <username> -p <password> get -r --downloadLocation data/cbioportal syn5521835
```

## Reporting bugs

Please double check and triple check your bug, and if you're super certain that we're not infallible, then [click here](https://github.com/Sage-Bionetworks/Genie/issues) to report the bug, along with a clear explanation of how we can find or reproduce it.

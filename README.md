# AACR Project GENIE

## Introduction

This repo documents code used to gather, QC, and standardize data uploaded by participating sites for AACR's Project GENIE (Genomics, Evidence, Neoplasia, Information, Exchange).

## Dependencies

These are tools and packages you will need, to be able to reproduce these results:

## Fetch Data

Using Synapse's [command-line client](http://python-docs.synapse.org/CommandLineClient.html), download all datasets into a subfolder named `data`:
```bash
synapse -u <username> -p <password> get -r --downloadLocation data syn5521835
```

## Reporting bugs

Please double check and triple check your bug, and if you're super certain that we're not infallible, then [click here](https://github.com/Sage-Bionetworks/Genie/issues) to report the bug, along with a clear explanation of how we can find or reproduce it.

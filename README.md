# py_shapeit

Snakemake recipe for running/interacting with SHAPEIT via a python environment

This module is immature and should be considered a rough beta. This was
developed for personal research use, I am placing online as I thought it
*may* be useful to others. No guarantees!

[SHAPEIT documentation](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

I use [scikit-allel](https://github.com/cggh/scikit-allel) package for working
with genetic variation data, also included is a script that converts SHAPEIT
output to an hdf5 file, which is a convenient way of handling large scale
genetic data.

## Recommended use

Clone this repository and rename as an analysis directory.

```
DIR=organismX_shapeit_date
git clone https://github.com/hardingnj/py_shapeit.git $DIR
```
This will create a directory containing the required files.

Edit files as needed. At the very least you will need to create a new `bam_locations.txt`, and edit the `config.yaml` file. Please feel free to push extensions of the `Snakefile` back to master.

## Overview of files

- *Snakefile*; the main file required by snakemake that encodes the pipeline.
- *config.yaml*; configuration, including filepaths etc.
- *submit.sh*; One possible example of how you may invoke snakemake.
- *bam_locations.txt*; describes where to find bam files for each sample in your vcf.


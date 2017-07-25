# snakeit

Snakemake recipe for running/interacting with SHAPEIT via a python environment

This workflow is immature and should be considered rough beta. 
This was developed for research use, I am placing online as I thought it *may* be useful to others. No guarantees!
[SHAPEIT documentation](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

I use the [scikit-allel](http://scikit-allel.readthedocs.io/en/latest/) module for working with genetic variation data, also included is a script that converts SHAPEIT output to an hdf5 file, 
which is a convenient way of handling large scale genetic data. For more about this see the above link.

## Recommended use

Follow the suggestions presented in the snakemake [deployment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) documentation.

ie:

```
# clone workflow into working directory
git clone https://github.com/hardingnj/snakeit.git path/to/workdir
cd path/to/workdir

# edit config and workflow as needed
vim config.yaml

# install dependencies into isolated environment
conda env create -n myworkflow --file environment.yaml

# activate environment
source activate myworkflow

# execute workflow
snakemake -n

```

Edit files as needed. At the very least you will need to create a new `bam_locations.txt`, and edit the `config.yaml` file. Please feel free to push extensions of the `Snakefile` back to master.

NOTE: To run the snakefile you will need `numpy` and `pandas` installed.

## Overview of files

- *Snakefile*; File that encodes the pipeline.
- *config.yaml*; Configuration, including filepaths etc.
- *submit.sh*; An example of how you may invoke snakemake.
- *bam_locations.txt*; describes where to find bam files for each sample in your vcf.
- *environment.yaml*; describes conda environment necessary to run tools.

## Overview of pipeline

The pipeline takes vcf files as an input, assuming one vcf per chromosome/contig.

Split these large vcfs into manageable chunks, as the extract PIRs step, and the shapeit phasing step are resource intensive.

Run extract PIRs on each chunk

Run shapeit on each chunk

Run ligate haplotypes on chunks to create a single phased output file

Run `shapeit_2_hdf5.py` to create a much easier to work with hdf5 file.

~ @hardingnj

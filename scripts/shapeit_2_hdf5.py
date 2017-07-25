#! /usr/bin/python

__author__ = 'Nicholas Harding'

import numpy as np
import h5py
import os.path
import gzip
from itertools import islice

if os.path.isfile(snakemake.output.hdf5):
    raise FileExistsError("outfile already exists")
else:
    h5file = h5py.File(snakemake.output.hdf5, mode="w")

chrom = h5file.create_group(snakemake.wildcards.chrom)

# Create the groups
grp_calldata = chrom.create_group("calldata")
grp_variants = chrom.create_group("variants")

# read samples file
fh_samples = gzip.open(snakemake.input.samples, 'rb')
samples_header = fh_samples.readline()
samples_desc = fh_samples.readline()
sample_info = [s.decode() for s in fh_samples.readlines()]

sample_names = np.array([s.rstrip().split(' ')
                         for s in sample_info], dtype="|S8")[:, 1]
n_sam = len(sample_names)

# count lines
number_sites = snakemake.params.max_sites
print("Max sites set at {0} snps.".format(number_sites))

# create objects
samples = chrom.create_dataset('samples', data=sample_names)

position = grp_variants.create_dataset('POS', (0, ),
                                       maxshape=(number_sites, ),
                                       dtype="int",
                                       compression="gzip",
                                       compression_opts=1)


identify = grp_variants.create_dataset('ID', (0, ),
                                       maxshape=(number_sites, ),
                                       dtype="S8",
                                       compression="gzip",
                                       compression_opts=1)


reference = grp_variants.create_dataset('REF', (0, ),
                                        maxshape=(number_sites, ),
                                        dtype="S1",
                                        compression="gzip",
                                        compression_opts=1)

alternate = grp_variants.create_dataset('ALT', (0, ),
                                        maxshape=(number_sites, ),
                                        dtype="S1",
                                        compression="gzip",
                                        compression_opts=1)

genotypes = grp_calldata.create_dataset('genotype', (0, n_sam, 2),
                                        maxshape=(number_sites, n_sam, 2),
                                        dtype="int",
                                        compression="gzip",
                                        compression_opts=1)

fh_haplotypes = gzip.open(snakemake.input.haplotypes, 'rb')

n = 0
print("loading haplotypes...")
while True:
    print(n, "read...")
    chunk = list(islice(fh_haplotypes, snakemake.params.chunk_size))
    if not chunk:
        break
    print("chunk has", len(chunk), "lines")
    as_np = np.array([line.rstrip().split(b' ') for line in chunk])
    print("chunk read and converted to numpy:", as_np.shape)

    position.resize(as_np.shape[0] + position.shape[0], axis=0)
    position[n:] = as_np[:, 2].astype('int')

    identify.resize(as_np.shape[0] + identify.shape[0], axis=0)
    identify[n:] = as_np[:, 1]

    alternate.resize(as_np.shape[0] + alternate.shape[0], axis=0)
    alternate[n:] = as_np[:, 4]

    reference.resize(as_np.shape[0] + reference.shape[0], axis=0)
    reference[n:] = as_np[:, 3]

    genotypes.resize(as_np.shape[0] + genotypes.shape[0], axis=0)
    genotypes[n:] = as_np[:, 5:].astype('int').reshape((-1, n_sam, 2))

    n += len(chunk)

h5file.close()
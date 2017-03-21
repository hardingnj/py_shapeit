#! /usr/bin/python

__author__ = 'Nicholas Harding'

import argparse
import numpy as np
import gzip
from itertools import islice
import h5py
import os.path

parser = argparse.ArgumentParser(description='Tool to convert shapeIt output '
                                             'to anhima readable')

parser.add_argument('haplotypes', help='haplotypes file')

parser.add_argument('samples', help='samples file')

parser.add_argument('-c', '--chunksize', dest='chunk_size', action='store',
                    type=int, default=100000,
                    help='Number of lines to read at once from file')

parser.add_argument('-z', '--compression', dest='comp_level', action='store',
                    type=int, default=1, help='zlib compression level')

parser.add_argument('--chr', dest='chr', action='store',
                    help='contig label to use as key at root of hdf5')

parser.add_argument('-o', '--out', dest='outfile', action='store',
                    help='path to write hdf5 file')

args = parser.parse_args()

if os.path.isfile(args.outfile):
    raise FileExistsError("outfile already exists")
else:
    h5file = h5py.File(args.outfile, mode="w")

chrom = h5file.create_group(args.chr)

# Create the groups
grp_calldata = chrom.create_group("calldata")
grp_variants = chrom.create_group("variants")

# read samples file
fh_samples = gzip.open(args.samples, 'rb')
samples_header = fh_samples.readline()
samples_desc = fh_samples.readline()
sample_info = [s.decode() for s in fh_samples.readlines()]

sample_names = np.array([s.rstrip().split(' ')
                         for s in sample_info], dtype="|S8")[:, 1]
n_sam = len(sample_names)

# count lines
number_sites = sum(1 for line in gzip.open(args.haplotypes))
print("Haplotypes file contains {0} snps.".format(number_sites))

# create objects
#'filters = Filters(complevel=args.comp_level, complib='zlib')
samples = chrom.create_dataset('samples', data=sample_names)

position = grp_variants.create_dataset('POS', (0, ),
                                       maxshape=(number_sites, ),
                                       dtype="int")


identify = grp_variants.create_dataset('ID', (0, ),
                                       maxshape=(number_sites, ),
                                       dtype="S8")


reference = grp_variants.create_dataset('REF', (0, ),
                                        maxshape=(number_sites, ),
                                        dtype="S1")

alternate = grp_variants.create_dataset('ALT', (0, ),
                                        maxshape=(number_sites, ),
                                        dtype="S1")

genotypes = grp_calldata.create_dataset('genotype', (0, n_sam, 2),
                                        maxshape=(number_sites, n_sam, 2),
                                        dtype="int")

fh_haplotypes = gzip.open(args.haplotypes, 'rb')

n = 0
while True:
    chunk = list(islice(fh_haplotypes, args.chunk_size))
    if not chunk:
        break
    as_np = np.array([line.rstrip().split(b' ') for line in chunk])

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
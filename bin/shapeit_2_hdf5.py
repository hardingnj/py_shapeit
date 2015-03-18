#! /usr/bin/python

__author__ = 'Nicholas Harding'

import argparse
from tables import Filters, openFile, IntAtom, StringAtom
import numpy as np
import gzip
from itertools import islice
import os

parser = argparse.ArgumentParser(description='Tool to convert shapeIt output '
                                             'to anhima readable')

parser.add_argument('haplotypes', help='haplotypes file')

parser.add_argument('samples', help='samples file')

parser.add_argument('-c', '--chunksize', dest='chunk_size', action='store',
                    type=int, default=1e5,
                    help='Number of lines to read at once from file')

parser.add_argument('-z', '--compression', dest='comp_level', action='store',
                    type=int, default=1, help='zlib compression level')

parser.add_argument('--chr', dest='chr', action='store',
                    help='contig label to use as key at root of hdf5')

parser.add_argument('-o', '--out', dest='outfile', action='store',
                    help='path to write hdf5 file')

args = parser.parse_args()

if os.path.isfile(args.outfile):
    raise IOError(args.outfile + ' already exists. Exiting...')

h5file = openFile(args.outfile, mode="w")
root = h5file.root

# Create the groups
chrom = h5file.create_group(root, args.chr)
grp_calldata = h5file.create_group(chrom, "calldata")
grp_variants = h5file.create_group(chrom, "variants")

# read samples file
fh_samples = gzip.open(args.samples, 'rb')
samples_header = fh_samples.readline()
samples_desc = fh_samples.readline()
sample_info = fh_samples.readlines()
sample_names = np.array([s.rstrip().split(' ') for s in sample_info])[:, 1]

# count lines
number_sites = sum(1 for line in gzip.open(args.haplotypes))
print "Haplotypes file contains {0} snps.".format(number_sites)

# create objects
filters = Filters(complevel=args.comp_level, complib='zlib')
samples = h5file.create_array(chrom, 'samples', sample_names)

position = h5file.create_earray(grp_variants, name='POS',
                                atom=IntAtom(itemsize=4),
                                expectedrows=number_sites, shape=(0, ),
                                filters=filters)

identify = h5file.create_earray(grp_variants, name='ID',
                                atom=StringAtom(itemsize=8),
                                expectedrows=number_sites, shape=(0, ),
                                filters=filters)

reference = h5file.create_earray(grp_variants, name='REF',
                                 atom=StringAtom(itemsize=1),
                                 expectedrows=number_sites, shape=(0, ),
                                 filters=filters)

alternate = h5file.create_earray(grp_variants, name='ALT',
                                 atom=StringAtom(itemsize=1),
                                 expectedrows=number_sites, shape=(0, ),
                                 filters=filters)

genotypes = h5file.create_earray(grp_calldata, name='genotype',
                                 atom=IntAtom(itemsize=1),
                                 expectedrows=number_sites,
                                 shape=(0, sample_names.size, 2),
                                 filters=filters)

fh_haplotypes = gzip.open(args.haplotypes, 'rb')

while True:
    chunk = list(islice(fh_haplotypes, args.chunk_size))
    if not chunk:
        break
    as_np = np.array([line.rstrip().split(' ') for line in chunk])
    try:
        position.append(as_np[:, 2].astype('int'))
        identify.append(as_np[:, 1])
        reference.append(as_np[:, 3])
        alternate.append(as_np[:, 4])

        geno = as_np[:, 5:].astype('int').reshape((-1, sample_names.size, 2))
        assert geno.shape[1] == sample_names.size
        genotypes.append(geno)

    except AssertionError:
        print "Assertion error, shapes don't match:", \
            sample_names.size, as_np.shape
        exit(1)

h5file.close()
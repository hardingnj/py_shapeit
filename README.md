# py_shapeit

Collection of functions for running/interacting w/ SHAPEIT from a python environment

This module is immature and should be considered a rough beta. This was
developed for personal research use, I am placing online as I thought it
*may* be useful to others. No guarantees!

Currently, the module assumes that jobs will be submitted using [SGE]
(http://star.mit.edu/cluster/docs/latest/guides/sge.html) and
dependencies are encoded with the hold_jid argument to qsub. In future I hope
 to use this module as a platform for experimenting with [luigi]
(https://github.com/spotify/luigi).

[SHAPEIT documentation](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

[Example code]("coming soon")

I use the [anhima](https://github.com/alimanfoo/anhima) package for working
with genetic variation data, also included is a script that converts SHAPEIT
output to an hdf5 file, which is a convenient way of handling large scale
genetic data.

## Installation
(After installing dependencies, tables, pyyaml)

```
git clone https://github.com/hardingnj/py_shapeit.git
cd py_shapeit
python setup.py install
```
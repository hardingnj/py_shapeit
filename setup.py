from ast import literal_eval
from distutils.core import setup


def get_version(source='src/shapeit/__init__.py'):
    with open(source) as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.split('=')[-1].lstrip())
    raise ValueError("__version__ not found")


setup(
    name='py_shapeit',
    version=get_version(),
    author='Nicholas Harding',
    author_email='njh@well.ox.ac.uk',
    package_dir={'': 'src'},
    packages=['shapeit'],
    url='https://github.com/hardingnj/py_shapeit',
    license='MIT License',
    description='Python wrapper to SHAPEIT tool',
    long_description=open('README.md').read(),
    classifiers=['Intended Audience :: Developers',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python',
                 'Topic :: Software Development :: Libraries :: Python Modules'
                 ],
    scripts=['bin/shapeit_2_hdf5.py'])
